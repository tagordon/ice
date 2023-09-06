clear all
close all

% load calculated DoS (to be implemented here to change with gamma)
load spEvib_VII
load all_iceVII_literature_data_number.mat

H2O_M = 18.01528 *1e-3 ; %kg/mol 

P_max= 800e3;
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
Vb=7.0441e-04;
x0 =([14e3 , Kp , .8e-3]); % K0, Kp, V
ifit=1;

Evib=fnval(spEvib, [data(:,2) (V(:)/Vb)]'); 
Pthermal=1e-6*gamma*V.^-1.*(V/Vb).^q.*Evib';

if ifit
    fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % 
    bestx = fminsearch(fun,x0);
    K0 = bestx(1);
    Vo = bestx(3);
    Kp = bestx(2);
% 
else
    fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % Variable Kp
    bestx = fminsearch(fun,x0([1 3]));
    K0 = bestx(1);
    Vo = bestx(2);
    Kp = x0(2);
end

Km=[Kp K0];

% "automatic" V, P parameters space boundary for plot and calculation
V_bound = ([.99*min(V) 1.01*max(V)]);
P_bound = ([min(data(:,2))-0.6*min(data(:,2)) max(data(:,2))+0.6*max(data(:,2))]);
Vc=linspace(V_bound(1),V_bound(2),200); 


Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]); % Cold compression curve
Pdc=finite_strain(V.^-1,Vo^-1,Km,[1 1]); % residual data - Cold compression curve
Pd=Pdc+Pthermal; % ideal data

%residual
sp_fit=csapi(Vc,Pc/1e3);
dP_data = data(:,1)/1e3-fnval(sp_fit,V);
dP_cold = data(:,1)/1e3-Pthermal(:)/1e3-fnval(sp_fit,V);

rms_cold = sqrt(sum((dP_cold).^2)/length(dP_cold));
rms_data = sqrt(sum((dP_data).^2)/length(dP_data));


% create lots of PVT points from the EOS determined above
Vc=linspace(0.9*V_bound(1),1.1*V_bound(2),100); 
Tc=linspace(0,T_max,100);
Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]);
Evib=fnval(spEvib, {Tc,Vc(:)/Vb});
spCv=fnder(spEvib,[1 0]); % will need the description of Cv(V.T)

% here is the calculation of the grid of PVT points
Pmm=Pc(:)*ones(1,length(Tc)) +1e-6*gamma*(Vc(:).^-1.*(Vc(:)/Vb).^q)*ones(length(Tc),1)'.*Evib';
Tmm=ones(length(Vc),1)*Tc(:)';
Vmm=Vc(:)*ones(1,length(Tc));

% put a LBF through V as a function of P and T - many more data than
% control points - no regularization needed.

npc = 10;  % adjust these to change the numbers of control points
ntc = 20;  % adjust these to change the numbers of control points
Pb = ([0, P_max]) ; % min and max Pressures
Tb = ([0, T_max]) ; % min and max temperature
Xc = linspace(Pb(1), Pb(2), npc);
Yc = linspace(Tb(1), Tb(2), ntc);

options.Xc = {Xc,Yc};
% full co-location solutions - no need for regularization
options.nReg=[1 1];
options.ordr=[6 6];
dY= 1e-4*Vmm(:);  % 1000ppm volume uncertainty?
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2) & Tmm(:)>=Tb(1) & Tmm(:)<Tb(2));  % choose the pressure range to fit.
sp_V_fPT=spdft([Pmm(id) Tmm(id)],Vmm(id),dY(id),options);

% start with defining a grid for calculations in P and T

% We're resetting the first value of T to eps, then relocating 
% the point closest to (P_ref, T_ref) at (P_ref, T_ref)
PT={0.:500:P_max,0.:10:T_max};
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
% diverges in the denominator
Cp= Cv- 1e6*Tm.*fnval(fnder(sp_V_fPT,[0 1]),PT).^2./fnval(fnder(sp_V_fPT,[1 0]),PT);
% make it nonnegative again (why eps rather than 0 though?)
id=find(Cp==0);
Cp(id)=eps;
rho=Vm.^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Cv from DoS model.  Use LBF for V(P,T) to calculate Cp.  at the
% end of this section have the ingredients for solution for Gibbs energy

% Plot V vs. P and T (here V is the interpolated version constructed 
% from the V, T grid above (with P determined by MG EOS)
figure(1)
subplot(1, 2, 1)
therm_surf(PT, Vm)
hold on
% compare this to the input data before interpolation
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2));
plot3(Pmm(id)*1e-3, Tmm(id), Vmm(id), 'ko','MarkerFaceColor','k','MarkerSize',10)
title("Mie-Grunensen P, V, T with Interpolation")
% and the residuals
subplot(1, 2, 2)
resid = Vmm - fnval(sp_V_fPT, {Pc, Tc})
plot3(Pmm(id)*1e-3, Tmm(id), resid(id), 'ko','MarkerFaceColor','k','MarkerSize',10)
title("Volume residuals")

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
npc=10;  % adjust these to change the numbers of control points
ntc=10;  % adjust these to change the numbers of control points
Pb = [ 0.  P_max ] ; % min and max Pressures
Tb = ([ 0 , T_max ]) ; % min and max temperature
Xc=linspace(Pb(1),Pb(2),npc);
Yc=linspace(Tb(1),Tb(2),ntc);
Yc=[0:2:20 30:10:200 logspace(log10(200),log10(T_max),5)];
clear Options
clear Data

Options.PTmc={Xc,Yc};
% full co-location solutions - no need forregularization
Options.ordr=[6 6];
Options.nReg=[];
%  Options.lam=1e0*[5e1 1e4];
Options.weight=[5e6 1e0];
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
figure(2)
subplot(1, 3, 1)
therm_surf(PT, fnval(sp_G_fPT, PT));
title("calculated gibbs surface")
hold on
out = fnGval(sp_G_fPT, PT);

subplot(1, 3, 2)
id = find(Pmm>=Pb(1) & Pmm<Pb(2) & Tmm>=Tb(1) & Tmm<Tb(2)); 
therm_surf(PT, 1./out.rho)
hold on
plot3(Pmm(id)*1e-3, Tmm(id), Vmm(id), 'ko','MarkerFaceColor','k','MarkerSize',3)
title("Volume from gibbs surface")

% residuals between P,V,T input and volume derived from gibbs surface
subplot(1, 3, 3)
id = find(Pmm>=Pb(1) & Pmm<Pb(2) & Tmm>=Tb(1) & Tmm<Tb(2)); 
plot3(Pmm(id)*1e-3, Tmm(id), Vmm(id)-1./fnGval(sp_G_fPT, {Pc, Tc}).rho(id), 'ko','MarkerFaceColor','k','MarkerSize',3)
title("volume residuals")

% residuals between other gibbs-derived properties

%Cp
figure(4)
subplot(1, 2, 1)
therm_surf(PT, out.Cp);
hold on
[Pg, Tg] = ndgrid(PT{1}, PT{2});
plot3(Pg*1e-3, Tg, Cp, 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Gibbs-derived Cp")
subplot(1, 2, 2)
plot3(Pg*1e-3, Tg, Cp - out.Cp, 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Cp Residuals")

%Cv
figure(5)
subplot(1, 2, 1)
therm_surf(PT, out.Cv);
hold on
[Pg, Tg] = ndgrid(PT{1}, PT{2});
plot3(Pg*1e-3, Tg, Cv, 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Gibbs-derived Cv")
subplot(1, 2, 2)
plot3(Pg*1e-3, Tg, Cv - out.Cv, 'ko','MarkerFaceColor','k','MarkerSize',3);
title("Cv Residuals")

% Get French for comparison 
load french_gibbs.mat

npc=10
ntc=10

Pb = ([0, 600e3])
Tb = ([10, 2000])
Xc = linspace(Pb(1), Pb(2), npc);
Yc = linspace(Tb(1), Tb(2), ntc)
options.Xc = {Xc, Yc}
options.nReg=[1, 1];
options.ordr=[6, 6];
options.lam=[1, 1];
dY=1d-3*G

id=find(p(:)>=Pb(1) & p(:)<Pb(2) & T(:)>=Tb(1) & T(:)<Tb(2));
gibbs_interp = spdft([transpose(p(id)), transpose(T(id))],G(id),dY(id),options);

PT={0.:100:600e3,10.:10:2000};
G_interp = fnval(gibbs_interp, PT);
minG = min(G_interp, [], 'all');

figure(3)
therm_surf(PT, (G_interp - minG).*1e6)
hold on 
plot3(p*1e-3, T, (G - minG).*1e6, 'ko','MarkerFaceColor','k','MarkerSize',10)
title("Normalized Gibbs Surface from French Model")

out = fnGval(gibbs_interp, PT);