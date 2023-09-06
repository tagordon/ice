load iceVII_m.mat

iceVII_m.Lin_2004.Simon(2)


dat_A =[3.48, 508.12
5.01, 555.22
7.51, 593.76
9.49, 733.99
14.01, 750.04
14.92, 727.56
17.40, 815.34
20.02, 914.90
21.01, 929.88
22.99, 1006.96
26.02, 1124.71];



P_N =[11.8, 1920.1
18.8, 2113.0
27.5, 2374.0
39.2, 2578.3
50.9, 2703.0
63.7, 2827.8
79.5, 2941.1
93.5, 3031.8
112.8, 3179.2
133.8, 3337.9
156.6, 3473.8
187.0, 3643.8
212.1, 3791.1
233.1, 3881.7
272.8, 4074.2
315.4, 4255.4
368.0, 4413.8
427.0, 4606.1];

P_U =[10.6, 2113.1
19.3, 2419.5
28.7, 2703.3
41.5, 2952.9
56.1, 3191.2
71.3, 3406.7
88.3, 3656.3
104.6, 3860.5
126.8, 4075.9
156.6, 4325.4
188.1, 4563.4
214.4, 4733.4
245.4, 4914.7
283.9, 5130.0
325.4, 5424.7
360.4, 5628.7
396.6, 5832.6
443.3, 5968.3
];


Melt_Millot = [116.3, 3610.6
135.0, 4053.1
162.4, 4518.3
187.0, 4904.1
221.4, 5255.6
260.5, 5698.0
311.9, 6219.6
370.3, 6718.5
422.3, 7126.7];



T_Lin=linspace(350,850,40);
P_Lin=((T_Lin/iceVII_m.Lin_2004.Simon(2)).^iceVII_m.Lin_2004.Simon(4)-1)*iceVII_m.Lin_2004.Simon(3)+iceVII_m.Lin_2004.Simon(1);


% iceVII (Wagner)
par_VII =[0.173683e1,0.544606e-1,0.806106e-7,2216,355];
res=.01;
MS=0;
TmVII = par_VII(5)-MS:res:700;
PmVII=par_VII(4)*exp( par_VII(1)*(1-(TmVII/par_VII(5)).^(-1)) - par_VII(2)*(1-(TmVII/par_VII(5)).^(5)) + par_VII(3)*(1-(TmVII/par_VII(5)).^(22)) );

Frank_lower=([
600 7.19
783 23.37
860 32.02
902 36.59]);


%plot(iceVII_m.Bridgman_1937.data(:,2),iceVII_m.Bridgman_1937.data(:,1),'o')
plot(iceVII_m.Datchi_2000.data(:,2),iceVII_m.Datchi_2000.data(:,1),'ko','MarkerFaceColor','b','MarkerSize',6)
hold on
plot(iceVII_m.Journaux_2013.data(:,2),iceVII_m.Journaux_2013.data(:,1),'d')
plot(iceVII_m.Frank_2004.data(:,2),iceVII_m.Frank_2004.data(:,1),'rv')
%plot(Frank_lower(:,2),Frank_lower(:,1),'rv')
plot(iceVII_m.Mishima_1978.data(:,2),iceVII_m.Mishima_1978.data(:,1),'o')
plot(iceVII_m.Pistorius.data(:,2),iceVII_m.Pistorius.data(:,1),'o')
plot(iceVII_m.Shwagner_2004.data(:,2),iceVII_m.Shwagner_2004.data(:,1),'.','MarkerSize',9)
plot(iceVII_m.Shwagner_2008.data(:,2),iceVII_m.Shwagner_2008.data(:,1),'.','MarkerSize',9)
plot(iceVII_m.Shwagner_2004.data_s(:,2),iceVII_m.Shwagner_2004.data_s(:,1),'r.')
plot(iceVII_m.Goncharov_2004.data(:,2),iceVII_m.Goncharov_2004.data(:,1),'r+')
plot(iceVII_m.Kimura_2014.data(:,2),iceVII_m.Kimura_2014.data(:,1),'b+')
plot (dat_A(:,1),dat_A(:,2),'ks','MarkerFaceColor','k','MarkerSize',6)
plot(PmVII.*1e-3,TmVII,'k-');
plot(P_Lin,T_Lin,'k--');


dat_K=[8.95 568.69 4.75
11.48 651.13 5.32
15.01 649.28 5.49
20.26 737.69 6.40
30.76 1302.74 7.30
37.19 1530.40 7.86
41.08 1567.78 7.85
47.32 1730.73 8.49
49.61 1766.10 8.52
51.46 1864.20 8.68
53.56 1911.32 8.89
55.11 1938.83 8.72
56.28 1946.71 8.87];


dat_Q=[8.39, 660.82
8.78, 676.89
11.30, 756.07
12.49, 790.49
14.58, 852.46
15.12, 905.25
16.59, 930.49
16.59, 941.97
17.28, 978.69
27.00, 1171.48
36.71, 1309.18
44.74, 1491.64];

plot (dat_K(:,1),dat_K(:,2),'ko','MarkerFaceColor','r','MarkerSize',9)

plot (dat_Q(:,1),dat_Q(:,2),'ks','MarkerFaceColor','r','MarkerSize',7)

% plot(187.0, 4904.1,'ko','MarkerFaceColor','b','MarkerSize',9)
% plot(Melt_Millot(:,1),Melt_Millot(:,2),'-','LineWidth',5)
% plot(P_N(:,1),P_N(:,2),'-','LineWidth',5)
% plot(P_U(:,1),P_U(:,2),'-','LineWidth',5)
% plot(187.0, 4904.1,'ko','MarkerFaceColor','b','MarkerSize',9)


plot([14.76 33.96], [842.13 507.05],'b--')
%plot([25.76 14.89], [644.75 497.87],'b--')
plot([19.51 43.33], [890 1150],'b--')
plot([49.79 61.04], [1767.23 1582.55],'r--')

legend('LBF','Datchi 2000 (optical)','Journaux 2013(optical)','Frank 2004 (XRD, only 110)'...
    ,'Mishima 1978 (resitance)',' Pistorius 1963 (resitance)','Shwagner 2004 (laser spekle patterns)',...
    'Shwagner 2008 (laser spekle patterns)', 'Shwagner 2008 (Solid-solid?)',...
    'Goncharov 2004 (Raman)','Kimura 2014 (Raman)', 'Ahart 2014 (Brillouin)','Wagner 2011 (IAPWS)','Lin 2004 (Raman)',...
    'Kimura 2020','Queyroux 2020','Millot et al. 2018 (shock)','melting line Millot et al. 2018 (shock)','Neptune adiabat','Uranus adiabat','Location','northwest')


% 
% 25.76, 644.75
% 14.89, 497.87
% 19.51, 966.07
% 43.33, 1227.70

text(40,1375,'\alpha-SI / VII"','FontSize',13)
text(65,1800,'\beta-SI / XVIII','FontSize',13)
text(30,700,'VII''','FontSize',13)
text(13,710,'Plastic?','FontSize',10)
text(20,500,'VII','FontSize',13)
text(70,800,'X','FontSize',13)
xlabel('Pressure (GPa)')
ylabel('Temperature (K)')



Tbound=[350 2500];
Pbound=[1 100];
ylim(Tbound)
xlim(Pbound)

