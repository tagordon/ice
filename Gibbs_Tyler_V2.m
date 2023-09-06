clear all
close all

% load calculated DoS (to be implemented here to change with gamma)
load spEvib_VII
load all_iceVII_literature_data_number.mat

H2O_M = 18.01528 *1e-3 ; %kg/mol 

P_max= 600e3;
T_max= 1000;

%Simon curve for Melting comparision (Wagner et al. 2011) :
P_ref = 2200; %  VI-VII transition
T_ref = 300; % room T

id=find(All_data(:,7) <= 16) ; 

% GPa units in the first sections

%  P (MPa)    T (K)    volume (m^3/kg)
data=[1e3*All_data(id,1)  All_data(id,2) All_data(id,3)*.1e-5./H2O_M];
V=data(:,3);

% convert volumes to specific volumes

% Data uncertainties
dV = 4e-3*V; % 0.2% in m^3/kg
dP = 0.05*data(:,1)*1e-3; % in GPa  (20 MPa nominal ruby uncertainty)

% Mie Grueisen EOS fit parameters
gamma=0.7;
q=1;
Kp=5;
Kpp=0.0;
Vb=7.0441e-04;
%x0 =([14e3 , Kp , .8e-3]); % K0, Kp, V
x0 = ([25.8e3, 4.2, -8e-5, 0.63e-3]); % K0, Kp, Kpp, V
ifit=1;

Evib=fnval(spEvib, [data(:,2) (V(:)/Vb)]'); 
Pthermal=1e-6*gamma*V.^-1.*(V/Vb).^q.*Evib';

fit_opts.MaxFunEvals = 1e4;
fit_opts.MaxIter = 1e4;

if ifit
    fun = @(x)MGval(x,Evib,gamma,data,dP,Kp,q,Vb); % 
    bestx = fminsearch(fun,x0,fit_opts);
    K0 = bestx(1);
    Vo = bestx(4);
    Kp = bestx(2);
    Kpp = bestx(3);
    Km=[Kpp Kp K0];
    %Km = [Kp, K0];
% 
else
    fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % Variable Kp
    bestx = fminsearch(fun,x0([1 3]));
    K0 = bestx(1);
    Vo = bestx(2);
    Kp = x0(2);
    Km=[Kp K0];
end

% "automatic" V, P parameters space boundary for plot and calculation
V_bound = ([.99*min(V) 1.01*max(V)]);
P_bound = ([min(data(:,2))-0.6*min(data(:,2)) max(data(:,2))+0.6*max(data(:,2))]);
Vc=linspace(V_bound(1),V_bound(2),200); 


Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]); % Cold compression curve
Pdc=finite_strain(V.^-1,Vo^-1,Km,[1 1]); % residual data - Cold compression curve
Pd=Pdc+Pthermal; % ideal data

%residual
sp_fit = csapi(Vc,Pc/1e3);
dP_data = data(:,1)/1e3-fnval(sp_fit,V);
dP_cold = data(:,1)/1e3-Pthermal(:)/1e3-fnval(sp_fit,V);

rms_cold = sqrt(sum((dP_cold).^2)/length(dP_cold));
rms_data = sqrt(sum((dP_data).^2)/length(dP_data));


% create lots of PVT points from the EOS determined above
Vc=linspace(0.7*V_bound(1),0.9*V_bound(2),200); 
%Vc=logspace(log10(0.7*V_bound(1)),log10(0.9*V_bound(2)),500);
Tc=linspace(0,T_max,40);
Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]);
Evib=fnval(spEvib, {Tc,Vc(:)/Vb});
spCv=fnder(spEvib,[1 0]); % will need the description of Cv(V.T)

% here is the calculation of the grid of PVT points
Pmm=Pc(:)*ones(1,length(Tc)) +1e-6*gamma*(Vc(:).^-1.*(Vc(:)/Vb).^q)*ones(length(Tc),1)'.*Evib';
Tmm=ones(length(Vc),1)*Tc(:)';
Vmm=Vc(:)*ones(1,length(Tc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put a LBF through V as a function of P and T - many more data than
% control points - no regularization needed.

npc = 30;  % adjust these to change the numbers of control points
ntc = 10;  % adjust these to change the numbers of control points
Pb = ([1000, P_max]) ; % min and max Pressures
Tb = ([0, T_max]) ; % min and max temperature
Xc = linspace(0, Pb(2), npc);
%Xc = [logspace(log10(Pb(1)), log10(200e3), 10) 250e3:50e3:Pb(2)];
%Xc = [0:10e3:500e3 510e3:10e3:600e3];
Yc = linspace(Tb(1), Tb(2), ntc);
%Yc = [Yc, [1820, 1840, 1860, 1880, 1900, 1920, 1940, 1960, 1980, 2000]];

options.Xc = {Xc,Yc};
% full co-location solutions - no need for regularization
%options.nReg=[2 2];
options.ordr=[6 6];
%options.lam=[1e2 1e2];
dY= 1e-6*Vmm(:);  % 1000ppm volume uncertainty?
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2) & Tmm(:)>=Tb(1) & Tmm(:)<Tb(2));  % choose the pressure range to fit.
sp_V_fPT=spdft([Pmm(id) Tmm(id)],Vmm(id),dY(id),options);

% start with defining a grid for calculations in P and T

% We're resetting the first value of T to eps, then relocating 
% the point closest to (P_ref, T_ref) at (P_ref, T_ref)
PT={1000.:100:P_max,0.:10:T_max};
PT{2}(1) = eps;
[M,id_Pref] = min((PT{1}-P_ref).^2);
PT{1}(id_Pref) = P_ref;
[M,id_Tref] = min((PT{2}-T_ref).^2);
PT{2}(id_Tref) = T_ref;

% interpolate V vs. P and T on fine grid 
Vm=fnval(sp_V_fPT,PT);
[Pm,Tm]=ndgrid(PT{1},PT{2});
[np,nt]=size(Pm);
% interpolate Cv (from Evib data) on grid 
Cv=fnval(spCv,[Tm(:) Vm(:)/Vo]');
% just force it to be nonnegative I guess 
id=find(Cv<0);
Cv(id)=0;
Cv=reshape(Cv,np,nt);
Cp = Cv - 1e6 * Tm .* fnval(fnder(sp_V_fPT,[0 1]),PT).^2 ./ fnval(fnder(sp_V_fPT,[1 0]),PT);
% make it nonnegative again (why eps rather than 0 though?)
id=find(Cp==0);
Cp(id)=eps;
rho=Vm.^-1;

% Plot V vs. P and T (here V is the interpolated version constructed 
% from the V, T grid above (with P determined by MG EOS)
figure(1)
subplot(1, 2, 1)
therm_surf(PT, Vm)
hold on
% compare this to the input data before interpolation
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2));
plot3(Pmm(id)*1e-3, Tmm(id), Vmm(id), 'ko','MarkerFaceColor','k','MarkerSize',3)
title("Mie-Grunensen P, V, T with Interpolation")
% and the residuals
subplot(1, 2, 2)
resid = Vmm - fnval(sp_V_fPT, {Pc, Tc});
plot3(Pmm(id)*1e-3, Tmm(id), resid(id), 'ko','MarkerFaceColor','k','MarkerSize',3)
title("Volume residuals")

%figure(1)
%plot3(Pm, Tm, Cp, 'ko','MarkerFaceColor','k','MarkerSize',3);
%title("Cp")

%figure(2)
%plot3(Pm, Tm, Cv, 'ko','MarkerFaceColor','k','MarkerSize',3);
%title("Cv")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LBF fit for Gibbs energy.  Trick is to center at ice VI triple point with
% ice V and water.

Cp_ref = Cp(id_Pref,:);

% Triple point between liquid ice III and ice V 
% get G0 from water EOS at triple point G0(Pref,Tref)
load('WaterEOS.mat'); % H2O LBF
load('IAPWS_sp_strct.mat');
sp_liq = G_H2O_2GPa_500K;
Ref_water = fnGval(sp_liq,{P_ref,T_ref});
S_water=Ref_water.S;
G_water=Ref_water.G;
dS=-3660.9-02.48-0.067943-19.384; %-3430 %-3660.9-02.48-0.067943 (IAPWS95)
S0 =S_water-dS; 
G_Pref= (S0.*(PT{2}-T_ref));

% set up spgft fit for gibbs energy 
npc=30;  % adjust these to change the numbers of control points
ntc=20;  % adjust these to change the numbers of control points
Pb = [ 1000.  P_max ] ; % min and max Pressures
Tb = ([ 1. , T_max ]) ; % min and max temperature
%Xc=linspace(220e3,Pb(2), npc);
Xc = [logspace(3,log10(100e3),15), 150e3:50e3:600e3];
%Xc=[linspace(Pb(1), 200e3, 10) Xc];
%Yc=linspace(100,Tb(2),ntc);
%Yc=[[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90] Yc];
Yc=[0:5:20 30:20:200 logspace(log10(200),log10(T_max),10)];
clear Options
clear Data

Options.PTmc={Xc,Yc};
% full co-location solutions - no need forregularization
Options.ordr=[6 6];
Options.nReg=[1 1];
%Options.lam=1e0*[1e8 1e2];
Options.weight=[1e2 1e0];
%Options.weight=1e2*[5e6 1e0];
Data.PTm=PT;

% corresponds to eq (2) in Journaux 2020 
G_ref = G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},Cp_ref./PT{2});
sp_G_ref = spaps(PT{2},G_ref,1);
deltaG_P_ref_T_ref = G_water - fnval(sp_G_ref,T_ref); 
Data.G = G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},Cp_ref./PT{2}) + deltaG_P_ref_T_ref;
 
Data.rho=rho;
Data.Cp=Cp; 
Data.MW = H2O_M;
Data.P_ref=P_ref;

% fit for the gibbs surface
sp_G_fPT=spgft_P(Data,Options);

% plot the gibbs surface 
figure(3)
%subplot(1, 3, 1)
therm_surf(PT, fnval(sp_G_fPT, PT));
title("calculated gibbs surface")
hold on
out = fnGval(sp_G_fPT, PT);

ds=1;

%Cp
figure(4)
subplot(1, 2, 1)
therm_surf(PT, out.Cp);
hold on
[Pg, Tg] = ndgrid(PT{1}, PT{2});
plot3(Pg(ds:ds:end)*1e-3, Tg(ds:ds:end), Cp(ds:ds:end), 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Gibbs-derived Cp with input")
subplot(1, 2, 2)
plot3(Pg(ds:ds:end)*1e-3, Tg(ds:ds:end), Cp(ds:ds:end) - out.Cp(ds:ds:end), 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Cp Residuals")

%Cv
figure(5)
subplot(1, 2, 1)
therm_surf(PT, out.Cv);
hold on
[Pg, Tg] = ndgrid(PT{1}, PT{2});
plot3(Pg(ds:ds:end)*1e-3, Tg(ds:ds:end), Cv(ds:ds:end), 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Gibbs-derived Cv with input")
subplot(1, 2, 2)
plot3(Pg(ds:ds:end)*1e-3, Tg(ds:ds:end), Cv(ds:ds:end) - out.Cv(ds:ds:end), 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Cv Residuals")