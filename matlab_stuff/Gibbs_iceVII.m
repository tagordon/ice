close all
clear all




%%
clear all
close all
% load calculated DoS (to be implemented here to change with gamma)
load spEvib_VII
load all_iceVII_literature_data_number.mat

H2O_M = 18.01528 *1e-3 ; %kg/mol 

%P_max=100000;
%T_max=1000;
P_max= 500e3;
T_max= 2200;

%Simon curve for Melting comparision (Wagner et al. 2011) :
P_ref = 2200; %  VI-VII transition
T_ref = 300; % room T

id=find(All_data(:,7) <= 16) ; 


% GPa units in the first sections

%  P (MPa)    T (K)    volume (m^3/kg)
data=[1e3*All_data(id,1)  All_data(id,2) All_data(id,3)*.1e-5./H2O_M];

V=data(:,3);
%

% convert volumes to specific volumes

% Data uncertainties
dV = 4e-3*V; % 0.2% in m^3/kg
dP = 0.05*data(:,1)*1e-3; % in GPa  (20 MPa nominal ruby uncertainty)

% Mie Grueisen EOS fit parameters
gamma=0.7;
q=1;
Kp=5;
Vb=7.0441e-04
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

figure(1)
%subplot(2,1,1)
errorbar(V(:),data(:,1)/1e3-Pthermal(:)/1e3,dP,dP,dV,dV,'kd','markerfacecolor','r')
hold on
errorbar(V,data(:,1)/1e3,dP,dP,dV,dV,'ko','markerfacecolor','w')
plot(V(:),data(:,1)/1e3,'ko','markerfacecolor','w','MarkerSize',6)
plot(V(:),data(:,1)/1e3-Pthermal(:)/1e3,'kd','markerfacecolor','r','MarkerSize',6)

plot(Vc,Pc/1e3,'k--','linewidth',2)
%plot(V,Pd,'ko','markersize',10)
hold off
xlabel('Volume (m^3/kg)')
ylabel('Pressure (GPa)')
title(sprintf('K0 = %0.1f GPa Kp = %0.1f, rho_0 = %0.0f kg/m^3 gamma = %0.1f',K0/1e3, Kp, Vo^-1, gamma))
%axis([3.8e-5 5.5e-5  -5 50])

%residual
sp_fit=csapi(Vc,Pc/1e3);
dP_data = data(:,1)/1e3-fnval(sp_fit,V);
dP_cold = data(:,1)/1e3-Pthermal(:)/1e3-fnval(sp_fit,V);

rms_cold = sqrt(sum((dP_cold).^2)/length(dP_cold));
rms_data = sqrt(sum((dP_data).^2)/length(dP_data));


%subplot(2,1,2)
plot(V,dP_data,'ko','markerfacecolor','w','MarkerSize',9)
hold on
plot(V,dP_cold,'kd','markerfacecolor','r')
plot([2e-4 8e-4],[0 0],'k-')

title(sprintf('rms data = %0.1f GPa ; rms (cold) = %0.1f GPa',rms_data,rms_cold))

%saveas(gcf,'EoSfit','jpg')
%saveas(gcf,'EoSfit','fig')

%%
% create lots of PVT points from the EOS determined above
%Vc=linspace(0.9*V_bound(1),1.1*V_bound(2),200); 
%Tc=linspace(0,T_max,100);
Vc=linspace(0.9*V_bound(1),1.1*V_bound(2),100); 
Tc=linspace(0,T_max,50);
Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]);
Evib=fnval(spEvib, {Tc,Vc(:)/Vb});
spCv=fnder(spEvib,[1 0]); % will need the description of Cv(V.T)

% here is the calculation of the grid of PVT points
Pmm=Pc(:)*ones(1,length(Tc)) +1e-6*gamma*(Vc(:).^-1.*(Vc(:)/Vb).^q)*ones(length(Tc),1)'.*Evib';
Tmm=ones(length(Vc),1)*Tc(:)';
Vmm=Vc(:)*ones(1,length(Tc));

%
% put a LBF through V as a function of P and T - many more data than
% control points - no regularization needed.

npc=20;  % adjust these to change the numbers of control points
ntc=20;  % adjust these to change the numbers of control points
Pb = ([ 0 , P_max ]) ; % min and max Pressures
Tb = ([ 0 , T_max ]) ; % min and max temperature
Xc=linspace(Pb(1),Pb(2),npc);
Yc=linspace(Tb(1),Tb(2),ntc);
clear options
options.Xc={Xc,Yc};
% full co-location solutions - no need for regularization
options.nReg=[1 1];
options.ordr=[6 6];
dY= 1e-4*Vmm(:);  % 1000ppm volume uncertainty?
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2));  % choose the pressure range to fit.
sp_V_fPT=spdft([Pmm(id) Tmm(id)],Vmm(id),dY(id),options);


%%
% calculate Cv from DoS model.  Use LBF for V(P,T) to calculate Cp.  at the
% end ofthis section have theingredients for solution forGibbs energy

% start with defining a grid for calculations in P and T
% This was a very dense grid (commented out below), which resulted 
% in a 32GB array. Trying a courser grid.
%PT={0.:5:P_max,0:.5:T_max};
PT={0.:500:P_max,0:10:T_max};
PT{2}(1)=eps;
 [M,id_Pref] = min((PT{1}-P_ref).^2);
 PT{1}(id_Pref)=P_ref;
 [M,id_Tref] = min((PT{2}-T_ref).^2);
 PT{2}(id_Tref)=T_ref;

Vm=fnval(sp_V_fPT,PT);
[Pm,Tm]=ndgrid(PT{1},PT{2});
[np,nt]=size(Pm);
Cv=fnval(spCv,[Tm(:) Vm(:)/Vo]');
id=find(Cv<0);
Cv(id)=0;
Cv=reshape(Cv,np,nt);
Cp= Cv- 1e6*Tm.*fnval(fnder(sp_V_fPT,[0 1]),PT).^2./fnval(fnder(sp_V_fPT,[1 0]),PT);
id=find(Cp==0);
Cp(id)=eps;
rho=fnval(sp_V_fPT,PT).^-1;

% P_ref profile in Cv
%   V_ref = fnval(sp_V_fPT,{P_ref,T_ref});
%   Cv_ref = fnval(spCv,{PT{2},V_ref/Vo});
%   Cp_ref = Cv_ref' - 1e6.* PT{2}.* fnval(fnder(sp_V_fPT,[0 1]),{P_ref,PT{2}}).^2./fnval(fnder(sp_V_fPT,[1 0]),{P_ref,PT{2}});

Cp_ref = Cp(id_Pref,:);

%
% LBF fit for Gibbs energy.  Trick is to center at ice VI triple point with
% ice V and water.

% Triple point between liquid ice III and ice V 

% get G0 from water EOS at triple point G0(Pref,Tref)
%Ref_water=IAPWS95({P_ref,T_ref},'P');

load('WaterEOS.mat'); % H2O LBF
load('IAPWS_sp_strct.mat');
sp_liq = G_H2O_2GPa_500K;
%sp_liq = sp_IAPWS;

Ref_water = fnGval(sp_liq,{P_ref,T_ref});



 
S_water=Ref_water.S;
G_water=Ref_water.G;

 
dS=-3660.9-02.48-0.067943-19.384; %-3430 %-3660.9-02.48-0.067943 (IAPWS95)
S0 =S_water-dS; 
G_Pref= (S0.*(PT{2}-T_ref));

 
npc=10;  % adjust these to chnage the numbers of control points
ntc=30;  % adjust these to chnage the numbers of control points
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

 
G_ref = G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},Cp_ref./PT{2});
sp_G_ref = spaps(PT{2},G_ref,1);
deltaG_P_ref_T_ref = G_water - fnval(sp_G_ref,T_ref); 
Data.G=  G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},Cp_ref./PT{2}) + deltaG_P_ref_T_ref;
 
Data.rho=rho;
Data.Cp=Cp; 
Data.MW = H2O_M;
Data.P_ref=P_ref;

 
sp_G_fPT=spgft_P(Data,Options);


%
figure(5)
out=fnGval(sp_G_fPT,PT);
subplot(211)
therm_surf(PT,1e6*(out.rho-Data.rho)/1e3,[],'\Delta rho (ppm)') %,[0 1000 240 400 -1e6 1e6])
subplot(212)
therm_surf({PT{1}(1:end),PT{2}},1e6*(out.Cp(1:end,:)-Data.Cp(1:end,:))/4e3)
zlabel('\Delta Cp (ppm)')



%%

PTplt={0:10:P_max,0:5:T_max};
mask_plt=ones(length(PTplt{1}),length(PTplt{2}));
out_plt=fnGval(sp_G_fPT,PTplt);


figure(2)
therm_surf(PTplt,out_plt.Kp,[],'Kp')
zlim([0 10])
view([33 -6])


iceIII_stability=([
  213.79484738501787 238.77551020408163 
   211.79950530816222 243.14868804664724
   209.83603770537775 248.2507288629738
   209.96353560166244 251.1661807580175
   238.44019073685286 252.332361516035
  269.0970598985118 253.3527696793003 
   304.13985669236445 254.66472303207001
   343.5558313287829 255.97667638483966
   345.58304787970974 252.332361516035
   343.2243367984429 248.39650145772595
  310.32350466217304 246.064139941691 
   273.0303699988949 243.29446064139944
   242.30975188909406 240.81632653061226
   213.79484738501787 238.77551020408163
]);
out_stb=fnGval(sp_G_fPT,[iceIII_stability]);

figure(3)
%
subplot(221)
therm_surf(PTplt,out_plt.rho,mask_plt,'Density')
shading flat
hold on
plot3(data(:,1),data(:,2),data(:,3).^-1,'ko','MarkerFaceColor','k','MarkerSize',10)
plot3(iceIII_stability(:,1),iceIII_stability(:,2), out_stb.rho,'r:','linewidth',2)
hold off
xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
zlabel('\rho (kg/m^{3})')
title('Ice III Density')
%
subplot(222)
%therm_surf(PTplt,out_plt.Cp,mask_plt,'Specific Heat')
%shading flat
therm_surf(PTplt,out_plt.Cv,mask_plt,'Specific Heat')
hold on
plot3(iceIII_stability(:,1),iceIII_stability(:,2), out_stb.Cv,'r:','linewidth',2)
hold off
xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
zlabel('Cv (J/kg/K)')
title('Ice III Specific Heat at constant volume')
%
subplot(223)
therm_surf(PTplt,out_plt.alpha,mask_plt,'Thermal Expansivity')%,[0 P_max 5 T_max 0 6e-4])
shading flat
hold on
plot3(iceIII_stability(:,1),iceIII_stability(:,2), out_stb.alpha,'r:','linewidth',2)
hold off
xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
zlabel('\alpha (K^{-1})')
title('Ice III Thermal Expansivity')
%
subplot(224)
therm_surf(PTplt,1e-3*out_plt.Kt,mask_plt,'Kt (GPa)')%,[0 P_max 5 T_max 0 6e-4])
shading flat
hold on
plot3(iceIII_stability(:,1),iceIII_stability(:,2), 1e-3*out_stb.Kt,'r:','linewidth',2)
hold off
xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
zlabel('Kt (GPa)')
title('Ice III Isothermal bulk modulus')

% therm_surf(PTplt,out_plt.vel,mask_plt,'Sound Speed')%,[0 P_max 0 T_max 2500 5000])
% shading flat

%saveas(gcf,'all_LBF','jpg')
%saveas(gcf,'all_LBF','fig')

%%
%
% now calculate melting phase boundary
out=fnGval(sp_G_fPT,PT);
outPTref=fnGval(sp_G_fPT,{P_ref,T_ref});
deltaS=S_water-outPTref.S;
mu_ice = H2O_M*(out.G);


% Melting data
%Data_phase_diagram_ices_Bridgman;
%Data_phase_diagram_ices_Grasset;

% Liquid water Gibbs values at same PT region for  form IAPWS
%Molal = 0;
%load('IAPWS_sp_strct.mat')
%water_results_1 = fnGval(sp_liq,PT);
%muw_water_1 = H2O_M*water_results_1.G;
%id=find(muw_water_1 == 0);
%muw_water_1(id) = NaN;
% extract PT coordinate of intersection of G_water and Gf :
%clear Melting_curve
%C = contours(PT{1},PT{2}, muw_water_1'-mu_ice', [0 0]);
%Melting_curve(:,1) = C(1,2:end);
%Melting_curve(:,2) = C(2,2:end);
%id=find(Melting_curve(:,1)>200);
%clear C
%figure(4)
%subplot(322)
%plot(P_melt,T_melt,'k-','linewidth',2)
%WaterPhaseDiagram_Wagner_2011
%hold on
%plot(Melting_curve(id,1),Melting_curve(id,2),'r--','linewidth',2)
%plot([P_ref,P_ref],[T_ref,220],'k-','linewidth',1)
%plot(P_ref,T_ref,'ko','MarkerSize',10)
%plot(350.1e-3,256.164,'ko','MarkerSize',10)
%plot([350.1,350.1],[256.164,220],'k-','linewidth',1)
%plot(TP_III_Bridgman(:,2)*1e3,TP_III_Bridgman(:,1),'o')
%plot(TP_V_Bridgman(:,2)*1e3,TP_V_Bridgman(:,1),'o')
%plot(TP_Ih_III_Bridgman(:,2)*1e3,TP_Ih_III_Bridgman(:,1),'o')
%plot(TP_III_V_Bridgman(:,2)*1e3,TP_III_V_Bridgman(:,1),'o')
%plot(TP_Ih_Bridgman(:,2)*1e3,TP_Ih_Bridgman(:,1),'o')

%plot(TP_III_Grasset(:,2),TP_III_Grasset(:,1),'*')
%plot(TP_Ih_Grasset(:,2),TP_Ih_Grasset(:,1),'*')
%plot(TP_V_Grasset(:,2),TP_V_Grasset(:,1),'*')

%xlabel('Pressure (MPa)')
%ylabel('Temperature (K)')
%title('ice III predicted melting curve (--) and pure H2O Simon melting curve (-) ')


%hold off 
%axis ([50 500 248 270])
%axis ([347 353 255.7 256.3])

%saveas(gcf,'Melting','jpg')
%saveas(gcf,'Melting','fig')

%%




%%
PTKs=[220 273.15-20 9.288 %"Tulk et al. (1997a)"
    250 273.15-20 9.7142
    280 273.15-20 9.781
    300 273.15-20 9.9581
    220 273.15-20 9.27 %"Tulk et al. (1997b)"
    276 273.15-27.2 9.94 
    210 273.15-25 8.5 %"Shaw et al. (1986)"
    340 273.15-25 9.3 %"Shaw et al. (1986)"
    276 273.15-27.2 9.87 % Gagnon et al. (1990)
    ];

out=fnGval(sp_G_fPT,PTKs(:,1:2));
disp("       P (MPa)     T (K)         Ks_ref        LBF Ks      residual (%)  ")
disp("      ----------------------------------------------")
disp([PTKs(:,1)  PTKs(:,2)    PTKs(:,3) out.Ks./1e3 (PTKs(:,3)-out.Ks./1e3)./PTKs(:,3)*100])

%% S0 residual :
format shortG
out_VI=fnGval(sp_G_fPT,{0.01,1});
S0_liqWa = -3516.467;  %J/K/kg
S0_theo = 189.9855597; %J/K/kg
S0_cal = out_VI.S - S0_liqWa;
S0_resi = out_VI.S - S0_liqWa - S0_theo ;
disp("S0 theoretical  S0 calculated   Residual")
disp([S0_theo    S0_cal               S0_resi])


%%
%save iceIII_IAPWS95_sp_G_fPT sp_G_fPT
%save iceIII_BolEoS_sp_G_fPT sp_G_fPT