tic

close all
clear all

% load calculated DoS (to be implemented here to change with gamma)
%load spEvib_VII_no_anarmonic.mat
load spEvib_VII.mat

flg_CCC =1; % 0: BM fit; 1: Helmotz LBF fit

ifit=1; % BM2(0) or BM3 (1)

H2O_M = 18.01528 *1e-3 ; %kg/mol 

P_max= 500e3;
T_max= 1000;

%triple point (Wagner et al. 2011) :
P_ref = 2.2e3;  %MPa
T_ref = 353.5; %K

Vo = 12.7218*1e-6/H2O_M

load all_iceVII_literature_data_number.mat

% convert volumes to specific volumes
id = find(All_data(:,7) <= 15);
data= [ 1e3*All_data(id,1) All_data(id,2)  All_data(id,3)*1e-6/H2O_M];
data = sortrows(data,3);
V = data(:,3);
% Data uncertainties
dV = 4e-3*V; % 0.2% in m^3/kg
dP = 0.00005*data(:,1); % in GPa  (20 MPa nominal ruby uncertainty)

% Mie Grueisen EOS fit parameters
gamma=1;
q=1;
Kp=6.5;
Vb=6.7069e-04;
x0 =([15e3 , Kp , 6.7069e-04]); % K0, Kp, V

Evib=fnval(spEvib, [data(:,2) (V(:)/Vb)]'); 
Pthermal=1e-6*gamma*V.^-1.*(V/Vb).^q.*Evib'; %MPa

V_bound = ([min(V)*0.95 max(V)*1.05]);
P_bound = ([min(data(:,2))-0.6*min(data(:,2)) max(data(:,2))+0.8*max(data(:,2))]);
Vc=linspace(V_bound(2),V_bound(1),200); 


if flg_CCC ==0
    if ifit
        fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % 
        bestx = fminsearch(fun,x0);
        K0 = bestx(1);
        Vo = bestx(3);
        Kp = bestx(2);
    else
        fun = @(x)MGval(x,Evib,gamma,data,Kp,q,Vb); % Variable Kp
        bestx = fminsearch(fun,x0([1 3]));
        K0 = bestx(1);
        Vo = bestx(2);
        Kp = x0(2);
    end 
    Km=[Kp K0];
    output
    Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]); % Cold compression curve
    Pdc=finite_strain(V.^-1,Vo^-1,Km,[1 1]); % residual data - Cold compression curve
else
    data_fit.PV=([data(:,1)-Pthermal data(:,3)]);
    % Sort data by increasing volume
    data_fit.PV = sortrows(data_fit.PV,2);
    % Define Vo
    data_fit.Vo = Vo ; % in cm^3/mol Klotz et al. (2017)
    options_ag.strainflg='log';
    options_ag.Vflg=1;
    np=18;
    options_ag.knt=linspace(12.73,3.7,np)*1e-6/H2O_M;
    options_ag.Reg=[linspace(12.73,11,np) 10 9 6  5.5 5 4.8 4.5 4 3.7]*1e-6/H2O_M;
    options_ag.Reg=[linspace(12.73,10.5,np) linspace(10,3.7,np)]*1e-6/H2O_M;
    options_ag.Reg=linspace(12.73,3.7,np)*1e-6/H2O_M;
    options_ag.drv=4;
    options_ag.lam=5e4;
    options_ag.k=6;
    %fit:
    sp_ag=sp_F_fit(data_fit,options_ag);
    % output
    out=fn_F_val(sp_ag,Vc.^-1); 
    Pc=out.P;
    out_c=fn_F_val(sp_ag,fliplr(V'.^-1));
    Pdc=fliplr(out_c.P);
end


% "automatic" V, P parameters space boundary for plot and calculation

Pd=Pdc+Pthermal; % ideal data

figure(1)
subplot(211)
errorbar(data(:,1)/1e3-Pthermal(:)/1e3,V(:),dV,dV,dP,dP,'rd','markerfacecolor','r')
hold on
errorbar(data(:,1)/1e3,V,dV,dV,dP,dP,'ko')
hold off

hold on
plot(Pc/1e3,Vc,'k--','linewidth',2)
%plot(V,data,'ko','markersize',10)
hold off
ylabel('Volume (m^3/kg)')
xlabel('Pressure (GPa)')
%title(sprintf('K0 = %0.1f GPa Kp = %0.1f, rho_0 = %0.0f kg/m^3 gamma = %0.1f',K0/1e3, Kp, Vo^-1, gamma))
%axis([V_bound(1) 7.8e-4 -0.30 3])

subplot(223)
plot(out.P*1e-3,out.K*1e-3,'k-','linewidth',1.5)
ylabel('Bulk Modulus (GPa)')
xlabel('Pressure (GPa)')

subplot(224)
plot(out.P*1e-3,out.Kp,'k-','linewidth',1.5)
ylabel('Kp')
xlabel('Pressure (GPa)')

saveas(gcf,'EoSfit','jpg')
saveas(gcf,'EoSfit','fig')

%%
% create lots of PVT points from the EOS determined above
Vc=linspace(V_bound(2),V_bound(1),50); 
Tc=linspace(0,T_max,50);

if flg_CCC ==0
    Pc=finite_strain(Vc.^-1,Vo^-1,Km,[1 1]);
else
    out=fn_F_val(sp_ag,Vc.^-1); 
    Pc=out.P;
end

figure
plot(Pc,Vc)
hold on
plot(data(:,1),data(:,3),'k.')
xlabel('Pressure (GPa)')
ylabel('Temperature (K)')



 
Evib=fnval(spEvib, {Tc,Vc(:)/Vb});
spCv=fnder(spEvib,[1 0]); % will need the description of Cv(V.T)

% here is the calculation of the grid of PVT points
Pmm=Pc(:)*ones(1,length(Tc)) +1e-6*gamma*(Vc(:).^-1.*(Vc(:)/Vb).^q)*ones(length(Tc),1)'.*Evib';
Tmm=ones(length(Vc),1)*Tc(:)';
Vmm=Vc(:)*ones(1,length(Tc));


% put a LBF through V as a function of P and T - many more data than
% control points - no regularization needed.

npc=20;  % adjust these to change the numbers of control points
ntc=10;  % adjust these to change the numbers of control points
Pb = ([ -1000 , P_max ]) ; % min and max Pressures
Tb = ([ 0 , T_max ]) ; % min and max temperature
Xc=linspace(Pb(1),Pb(2),npc);
Xc=[Pb(1) linspace(0,50e3,10) linspace(60e3,Pb(2),10)];
Yc=linspace(Tb(1),Tb(2),ntc);
clear options
options.Xc={Xc,Yc};
% full co-location solutions - no need for regularization
options.nReg=[3 2];
options.lam=[1e8 1e3];
options.ordr=[5 4];
dY= 1e-3*Vmm(:);  % 10th percent volume uncertainty?
id=find(Pmm(:)>=Pb(1) & Pmm(:)<Pb(2));  % choose the pressure range to fit.
sp_V_fPT=spdft([Pmm(id) Tmm(id)],Vmm(id),dY(id),options);

%
figure(5)
fnplt(sp_V_fPT)
hold on
plot3(data(:,1),data(:,2),data(:,3),'ko')
xlabel('Pressure (GPa)')
ylabel('Temperature (K)')
zlabel('Volume')
hold off

figure(6)
fnplt(fnder(sp_V_fPT,[1 0]))
hold off
%
%Plot residual
V_data_residual = (data(:,3)-fnval(sp_V_fPT,[data(:,1) , data(:,2)]')')./data(:,3).*100;

figure(7)
subplot(211)
plot(data(:,1)*1e-3,V_data_residual,'ko')
xlabel('Pressure (GPa)')
ylabel('Volume residual (%)')
grid on
subplot(212)
plot(data(:,2),V_data_residual,'ko')
xlabel('Temperature (K)')
ylabel('Volume residual (%)')
ylim =([-10 10])
grid on
hold off
%%
toc 
%%
% calculate Cv from DoS model.  Use LBF for V(P,T) to calculate Cp.  at the
% end ofthis section have theingredients for solution forGibbs energy
tic
% start with defining a grid for calculations in P and T
PT={0.1:100:P_max,0.01:20:T_max};
%PT{2}(1)=eps;
[M,id_Pref] = min((PT{1}-P_ref).^2);
PT{1}(id_Pref)=P_ref;
[M,id_Tref] = min((PT{2}-T_ref).^2);
PT{2}(id_Tref)=T_ref;

Vm=fnval(sp_V_fPT,PT);
id=find(Vm>Vo);
Vm(id)=Vo;
[Pm,Tm]=ndgrid(PT{1},PT{2});
[np,nt]=size(Pm);
Cv=fnval(spCv,[Tm(:) Vm(:)/Vo]');
id=find(Cv<=0);
Cv(id)=eps;
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
%%
figure; 
surf(PT{2},PT{1},Cp,'linestyle','none')
toc
%%
% LBF fit for Gibbs energy.  Trick is to center at ice VI triple point with
% ice V and water.

% Triple point between liquid ice III and ice V 

% get G0 from water EOS at triple point G0(Pref,Tref)
load('WaterEOS.mat'); % H2O LBF
load('IAPWS_sp_strct.mat');
sp_liq = G_H2O_2GPa_500K;
sp_liq = sp_IAPWS;

Ref_water = fnGval(sp_liq,{P_ref,T_ref});

 
S_water=Ref_water.S;
G_water=Ref_water.G;

  
dS=-3529.5+434.3;
S0 =S_water-dS; 
G_Pref= (S0.*(PT{2}-T_ref));

 
npc=9;  % adjust these to chnage the numbers of control points
ntc=12;  % adjust these to chnage the numbers of control points
Pb = [ 0.  P_max ] ; % min and max Pressures
Tb = ([ 0.01 , T_max ]) ; % min and max temperature
Xc=linspace(Pb(1),Pb(2),npc);
Yc=linspace(Tb(1),Tb(2),ntc);
%Yc=[0:2:50 50:5:200 logspace(log10(200),log10(T_max),20)];

clear Options
clear Data
Options.PTmc={Xc,Yc};
% full co-location solutions - no need forregularization
Options.ordr=[8 8];
Options.nReg=[ ];
 Options.lam=1e0*[5e12 1e4];
 Options.weight=1e0*[5e7 1e0];
Data.PTm=PT;

CpoT=Cp_ref./PT{2};
CpoT(1)=eps;
 
G_ref = G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},(CpoT));
sp_G_ref = spaps(PT{2},G_ref,1);
deltaG_P_ref_T_ref = G_water - fnval(sp_G_ref,T_ref); 
Data.G=  G_Pref + cumtrapz(PT{2},Cp_ref)-PT{2}.*cumtrapz(PT{2},CpoT) + deltaG_P_ref_T_ref;


 
Data.rho=rho;
Data.Cp=Cp; 
Data.MW = H2O_M;
Data.P_ref=P_ref;

sp_G_fPT=spgft_P(Data,Options);

figure 
%

% G0=fnval(sp_G_fPT,{P_ref,T_ref}); 
% sp_G_fPT.coefs=sp_G_fPT.coefs-G0+G_water;
%
figure(5)
out=fnGval(sp_G_fPT,PT);
subplot(211)
therm_surf(PT,1e6*(out.rho-Data.rho)/1e3,[],'\Delta rho (ppm)') %,[0 1000 240 400 -1e6 1e6])
subplot(212)
therm_surf({PT{1}(1:end),PT{2}},1e2*(out.Cp(1:end,:)-Data.Cp(1:end,:))/4e3)
zlabel('\Delta Cp (%)')


%

PTplt={0:2000:P_max,0:10:T_max};
mask_plt=ones(length(PTplt{1}),length(PTplt{2}));
out_plt=fnGval(sp_G_fPT,PTplt);


figure(2)
therm_surf(PTplt,out_plt.Kp,[],'Kp')
zlim([0 10])
view([33 -6])

figure(1)
subplot(221)
therm_surf(PTplt,out_plt.rho,mask_plt,'Density')
shading flat
hold on
plot3(1e0*data(:,1),data(:,2),data(:,3).^-1,'ko','MarkerFaceColor','k','MarkerSize',10)
hold off
subplot(222)
therm_surf(PTplt,out_plt.Cp,mask_plt,'Specific Heat')
shading flat
hold on
therm_surf(PTplt,out_plt.Cv,mask_plt,'Specific Heat')
hold off
%
subplot(223)
therm_surf(PTplt,out_plt.alpha,mask_plt,'Thermal Expansivity')%,[0 P_max 5 T_max 0 6e-4])
shading flat
%
subplot(224)
therm_surf(PTplt,1e-3*out_plt.Kt,mask_plt,'Kt (GPa)')%,[0 P_max 5 T_max 0 6e-4])
shading flat

% therm_surf(PTplt,out_plt.vel,mask_plt,'Sound Speed')%,[0 P_max 0 T_max 2500 5000])
% shading flat

saveas(gcf,'all_LBF','jpg')
saveas(gcf,'all_LBF','fig')

%
% now calculate melting phase boundary
out=fnGval(sp_G_fPT,PT);
outPTref=fnGval(sp_G_fPT,{P_ref,T_ref});
deltaS=S_water-outPTref.S;
mu_ice = H2O_M*(out.G);

% Melting data
Data_phase_diagram_ices_Bridgman;
Data_phase_diagram_ices_Grasset;


% Liquid water Gibbs values at same PT region for  form IAPWS
Molal = 0;
water_results_1 = fnGval(sp_liq,PT);
muw_water_1 = H2O_M*water_results_1.G;
id=find(muw_water_1 == 0);
muw_water_1(id) = NaN;
% extract PT coordinate of intersection of G_water and Gf :
clear Melting_curve
C = contours(PT{1},PT{2}, muw_water_1'-mu_ice', [0 0]);
Melting_curve(:,1) = C(1,2:end);
Melting_curve(:,2) = C(2,2:end);
id=find(Melting_curve(:,1)>300);
clear C
figure(4)
%subplot(322)
%plot(P_melt,T_melt,'k-','linewidth',2)
WaterPhaseDiagram_Wagner_2011
hold on
plot(Melting_curve(id,1),Melting_curve(id,2),'r--','linewidth',2)
plot([P_ref,P_ref],[T_ref,220],'k-','linewidth',1)
plot(P_ref,T_ref,'ko','MarkerSize',10)
plot(2216,355,'ko','MarkerSize',10) % Wagner et al. 2011
plot(2170,354.8,'go','MarkerSize',10) % Datchi 2000
plot(2190,354.8,'bo','MarkerSize',10) % Bridgman 1937
plot(2160,352.2,'ro','MarkerSize',10) % Journaux 2013
plot(2160,352.2,'ro','MarkerSize',10) % Journaux 2013

plot([632.4,632.4],[273.31,220],'k-','linewidth',1)
plot(TP_V_Bridgman(:,2)*1e3,TP_V_Bridgman(:,1),'o')
plot(TP_VI_Bridgman(:,2)*1e3,TP_VI_Bridgman(:,1),'o')
plot(TP_III_V_Bridgman(:,2)*1e3,TP_III_V_Bridgman(:,1),'o')
plot(TP_V_VI_Bridgman(:,2)*1e3,TP_V_VI_Bridgman(:,1),'o')
plot(TP_VI_Grasset(:,2),TP_VI_Grasset(:,1),'*')
plot(TP_V_Grasset(:,2),TP_V_Grasset(:,1),'*')
xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
title('ice VI predicted melting curve (--) and pure H2O Simon melting curve (-) ')

hold off 
% txt=sprintf('Delta S = %0.0f J/kg/K',deltaS);
% text(500, 370,txt)
% t = text(1500,290,'VI');
% s = t.FontSize;
% t.FontSize = 15;
% t = text(800,350,'Liquid H_2O');
% s = t.FontSize;
% t.FontSize = 15;
% t = text(2250,300,'VII');
% s = t.FontSize;
% t.FontSize = 15;

%axis ([350 2400 250 380])
axis ([2050 2300 340 365])
saveas(gcf,'Melting','jpg')
saveas(gcf,'Melting','fig')


%% COmpute bulk modulus from elastic moduli

Shimizu6 = [
1.177264   33.136966    27.746686   14.491900   11.134021   6.097202     5.390280
1.219653   32.783505  27.746686    14.580265    11.575847    5.920471    5.478645
1.373796   33.755523  28.188513    15.022091  12.989691    6.273932   5.920471
1.396917   34.108984  28.630339    15.110457   12.547865    6.097202    5.655376
1.601156   34.815906   30.397644   15.905744    12.812960   6.450663    5.832106
1.778420   36.494845  32.783505    17.673049    14.315169   6.008837     6.539028
2.005780   38.969072   33.667158    18.468336    14.933726    6.804124     6.450663
2.109827   39.941090  34.197349    18.733432   14.315169   6.980854   6.627393
 ];
Shimizu6=Shimizu6(:,[1 2 5 4 3 6 7]) ; % entered in ?wrong? order, shifted to cyclic C11 C12 C13 C33 etc

 
tulk6=[
    .62  1342.7 25.9 14.8 12.1 25.5 6.5 10.1
    .65 1347.8 26.4 14.6 12.8 25.8 6.4 10.1
    .72 1353.5 26.8 14.5 12.8 26.2 6.3 10.4
    .77 1359.6 27.5 15.5 13.3 26.8 6.2 10.6
    .82 1362.9 27.9 15.9 13.4 27 6.4 10.7
    ];

 
Ks6=zeros(8,1);
mus6=Ks6;
P6=Shimizu6(:,1);
for i=1:8
    hs=HSBounds(Ci2Cij(Shimizu6(i,2:end),'tet6'));
    Ks6(i)=mean(hs(:,1));
    mus6(i)=mean(hs(:,2));
end

 
Kt6=zeros(5,1);
mut6=Kt6;
P6t=tulk6(:,1);
for i=1:5
    hs=HSBounds(Ci2Cij(tulk6(i,3:end),'tet6'));
    Kt6(i)=mean(hs(:,1));
    mut6(i)=mean(hs(:,2));
end


%%
%     620 273.15-20 17.206 % Tulk et al. (1997) calculated with the wrong formula
%     650 273.15-20 17.636 % Tulk et al. (1997)
%     720 273.15-20 17.75 % Tulk et al. (1997)
%     770 273.15-20 18.401 % Tulk et al. (1997)
%     820 273.15-20 18.631 % Tulk et al. (1997)
PTKs=[[P6.*1e3  ones(length(P6),1).*(300) Ks6];
[P6t.*1e3 ones(length(P6t),1).*(273.15-20) Kt6];
[    
    720 273.15-20 17.77 % Tulk et al. (1997) 
    720 273.15-35 17.93 % Tulk et al. (1997)
    
    620 273.15-25 14.9 % Shaw et al. (1986) 
    800 273.15-25 16 % Shaw et al. (1986) 
    777 273-35.4 18.14  % Gagnon et al. (1990) 
    1230 300 19.5 % Shimizu (1996) 
    1700 300 22.4 % Baer (1998) 
    ]];


out=fnGval(sp_G_fPT,PTKs(:,1:2));
disp("       P (MPa)     T (K)         Ks_ref        LBF Ks      residual (%)  ")
disp("      ----------------------------------------------")
disp([PTKs(:,1)  PTKs(:,2)    PTKs(:,3) out.Ks./1e3 (PTKs(:,3)-out.Ks./1e3)./PTKs(:,3)*100])


out2=fnGval(sp_G_fPT,{600:10:2400,300});

figure
plot(PTKs(:,1),PTKs(:,3),'s')
hold on
plot(PTKs(:,1),out.Ks./1e3,'d')
plot(600:10:2400,out2.Ks./1e3,'-')
hold off



%% S0 residual :
format shortG
out_VI=fnGval(sp_G_fPT,{0.01,1});
S0_liqWa = -3516.467;  %J/K/kg
S0_theo = 194.4739415; %J/K/kg
S0_cal = out_VI.S - S0_liqWa;
S0_resi = out_VI.S - S0_liqWa - S0_theo ;
disp("S0 theoretical  S0 calculated   Residual")
disp([S0_theo    S0_cal               S0_resi])

%%
save iceVI_IAPWS95_sp_G_fPT sp_G_fPT