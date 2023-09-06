clear all
close all

% load calculated DoS (to be implemented here to change with gamma)
load spEvib_VII
load all_iceVII_literature_data_number.mat
load Data_litterature.mat

P_max= 800e3;
T_max= 1000;

id=find(All_data(:,7) <= 16);
H2O_M = 18.01528 *1e-3 ; %kg/mol 

%  P (MPa)    T (K)    volume (m^3/kg)
data=[1e3*All_data(id,1)  All_data(id,2) All_data(id,3)*.1e-5./H2O_M];
V=data(:,3);
P=data(:,1);
T=data(:,2);

% Data uncertainties
dV = 4e-3*V; % 0.2% in m^3/kg
dP = 0.05*data(:,1)*1e-3; % in GPa  (20 MPa nominal ruby uncertainty)

% Mie Grueisen EOS fit parameters
gamma=0.7;
q=1;
Vo = 0.63e-3;
%Kp=3.8;
%Kpp=0.0;
Vb=7.0441e-04;
%x0 =([14e3 , Kp , .8e-3]); % K0, Kp, V
%x0 = ([19.2e3, Kp, -0.09e-3, 0.63e-3]); % K0, Kp, Kpp, V
%ifit=1;

%Evib=fnval(spEvib, [data(:,2) (V(:)/Vb)]'); 
%Pthermal=1e-6*gamma*V.^-1.*(V/Vb).^q.*Evib';

%fit_opts.MaxFunEvals = 1e4;
%fit_opts.MaxIter = 1e4;

%if ifit
%    fun = @(x)MGval(x,Evib,gamma,data,dP,Kp,q,Vb); % 
%    bestx = fminsearch(fun,x0,fit_opts);
%    K0 = bestx(1);
%    Vo = bestx(4)
%    Kp = bestx(2);
%    Kpp = bestx(3);
%    Km=[Kpp Kp K0];
%    %gamma = bestx(5);
%    %Km = [Kp, K0];
%% 
%else
%    fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % Variable Kp
%%    bestx = fminsearch(fun,x0([1 3]));
%    K0 = bestx(1);
%    Vo = bestx(2);
%    Kp = x0(2);
%    Km=[Kp K0];
%end

%Km = [-8e-5, 4.2, 25.8e3]
%Vo = 0.63e-3

%Km = [-0.09e-3, 3.8, 19.2e3];

% take data at around 300K and remove  Hemley (1987) and Loubeyre (1999)
% from the fit
id=find(Data.VTP.data(:,2) <= 310  &  Data.VTP.data(:,4) ~= 3 & Data.VTP.data(:,4) ~= 8);
datap.PV=[Data.VTP.data(id,[3,1])];

% Sort data by increasing volume
datap.PV = sortrows(datap.PV,2);
% Define Vo
datap.Vo = 12.7218 ; % in cm^3/mol Klotz et al. (2017)



%%
%%%%%%%%%%% 4th order global fit
clear options
options.strainflg='Eulerian';
options.Vflg=1; %input in volumes
options.knt=[12.75  3.7]; %knots in volumes
options.k=5; %specify spline order (degree is one less)
sp_ref=sp_F_fit(datap,options);

% Compute volumes associated with the bulk modulus data
Vcp=linspace(12,4,1000);
out=fn_F_val(sp_ref,Vcp.^-1);

%[V, I] = sort(V); 
%outm = fn_F_val(sp_ref, V.^-1);
%P = data(:,1);
%Pm = P(I);

% "automatic" V, P parameters space boundary for plot and calculation
V_bound = ([.99*min(V) 1.01*max(V)]);
P_bound = ([min(data(:,2))-0.6*min(data(:,2)) max(data(:,2))+0.6*max(data(:,2))]);
Vc=linspace(V_bound(1),V_bound(2),200); 

%[Pc, Kc, Gc, Lc, Fc, Uc]=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]); % Cold compression curve
%[Pdc, Kdc, Gdc, Ldc, Fdc, Udc]=finite_strain(V.^-1,Vo^-1,Km,[1 1]); % residual data - Cold compression curve
%Pd=Pdc+Pthermal; % ideal data

%plot(Vc,Pc/1e3,'k--','linewidth',2)

%if (length(x0)==4)
%    title(sprintf('K0 = %0.1f GPa Kp = %0.1f, Kpp= %0.1fe-6, rho_0 = %0.0f kg/m^3 gamma = %0.1f',K0/1e3, Kp, Kpp*1e6, Vo^-1, gamma))
%else
%    title(sprintf('K0 = %0.1f GPa Kp = %0.1f, rho_0 = %0.0f kg/m^3 gamma = %0.1f',K0/1e3, Kp, Vo^-1, gamma))
%end

Vs = [];
Ps = [];
Ts = [];
Pres = [];
id=find(All_data(:,7) <= 16);
U = unique(All_data(id,end));
for i = 1:length(U)
    refid = find(All_data(id,end)==U(i));
    [Vtmp, Is] = sort(V(refid));
    Ptmp = P(refid);
    Ttmp = T(refid);
    Vs = [Vs; Vtmp];
    Ps = [Ps; Ptmp(Is)];
    Ts = [Ts; Ttmp(Is)];
    outseg = fn_F_val(sp_ref,flip(Vtmp*18.01528e3).^-1);
    Pres = [Pres; flip(outseg.P)];
end
dP = 0.05*Ps*1e-3; % in GPa  (20 MPa nominal ruby uncertainty)
dV = 4e-3*Vs; % 0.2% in m^3/kg

Evib=fnval(spEvib, [Ts(:) (Vs(:)/Vb)]'); 
Ptherm=1e-6*gamma*Vs.^-1.*(Vs/Vb).^q.*Evib';

figure(1)
subplot(2,1,1)
errorbar(Vs,Ps/1e3-Ptherm/1e3,dP,dP,dV,dV,'kd','markerfacecolor','r')
hold on
errorbar(Vs, Ps/1e3,dP,dP,dV,dV,'ko','markerfacecolor','w')
plot(Vs,Ps/1e3,'ko','markerfacecolor','w','MarkerSize',6)
plot(Vs,Ps/1e3-Ptherm/1e3,'kd','markerfacecolor','r','MarkerSize',6)

plot(Vcp/18.01528e3, out.P, 'k--', 'linewidth',2)
%plot(V,Pd,'ko','markersize',10)
hold off
xlabel('Volume (m^3/kg)')
ylabel('Pressure (GPa)')


subplot(2, 1, 2)
errorbar(Vs,Ps/1e3 - Pres,dP,dP,dV,dV,'ko','markerfacecolor','w')
hold on
plot(Vs, Ps/1e3 - Ptherm/1e3 - Pres,'kd','markerfacecolor','r','MarkerSize',6)
xlabel('Volume (m^3/kg)')
ylabel('Pressure residuals (GPa)')

figure(2)
% plot all bulk modulus data
plot(Data.K.Asahara2010.K_comp(:, 1), Data.K.Asahara2010.K_comp(:, 2), 'k.')
hold on
plot(Data.K.Shimizu1995.P(:), Data.K.Shimizu1995.Bs(:), 'b.')
%plot(Data.K.Shimizu1995.C11(:, 1), Data.K.Shimizu1995.C11(:, 2), 'b.')
plot(Data.K.Ahart2011(:, 1), Data.K.Ahart2011(:, 2), 'r.')
plot(Data.K.Zang2019.All(:, 1), Data.K.Zang2019.All(:, 5), 'g.')
plot(out.P(:), out.K(:), 'k-')
xlabel("Pressure (Gpa)")
ylabel("K (Gpa)")
title("Bulk modulus data")

figure(3)
subplot(2, 1, 1)
load('dos_high_pressure_ice/dos.mat');
[Tc, Vc] = ndgrid(sort(dos.T_array), sort(dos.VovVo));
surf(Tc, Vc, fnval(spEvib, {sort(dos.T_array), sort(dos.VovVo)}))
shading flat
hold on
plot3(dos.T_array, dos.VovVo, dos.Evib_array, 'ko')
ylabel('V/Vo')
xlabel('Temperature (K)')
title("LBF fit to Evib with data")
subplot(2, 1, 2)
plot3(dos.T_array, dos.VovVo, dos.Evib_array - fnval(spEvib, [dos.T_array; dos.VovVo]), 'ko')
xlabel('Temperature (K)')
ylabel('V/Vo')
title("Residuals")