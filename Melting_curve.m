clear all
close all

load('gibbs_vii.mat');
%load('iceVII_m.mat');

P = linspace(1.8e9, 100e9, 1000);
T = linspace(20, 1800, 1000);
PT = {P./1e6, T};
ice_vii = fnGval(gibbs_interp, PT);
water = SeaFreeze(PT, 'water2');

% find G0 at triple point
%water_triple = SeaFreeze([2200, 353.5], 'water2');
%gibbs_triple = fnGval(gibbs_interp, {[2200], [353.5]});
%G0 = water_triple.G - gibbs_triple.G
%ice_vii.G = ice_vii.G + G0;

% correct with S0
%S0 = -3100;
%G_corrected = ice_vii.G - S0 * (T - 353.5);
%Sc = fnGval(gibbs_interp, {[2200], [20]}).S - S0;

melt_line = contour(P, T, water.G' - ice_vii.G', [0, 0]);

figure(1)
%plot(melt_line(1,1600:2900)/1e9, melt_line(2,1600:2900), 'g-', linewidth=3)
plot(melt_line(1,:)/1e9, melt_line(2,:), 'g-', linewidth=3)
hold on
plot_data_melt

%figure(2)
%therm_surf(PT, ice_vii.G);
%hold on
%therm_surf(PT, water.G);