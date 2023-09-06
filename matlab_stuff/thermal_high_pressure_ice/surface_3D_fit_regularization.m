%% Import data
% close all
% clear all

%addpath('./Thermo_LBF_Tools/');

%% Evaluating spline for ion-thermal internal energy with regularization
%compute the internal energy from vibrational DOS data (comment if this has
%already been run once)
% high_pressure_ice_DOS;

% Check input data
n = 50;
X_data = [dos.T_array zeros(1,n) ones(1,n)*150.0];
Y_data = [dos.VovVo_array  linspace(min(dos.VovVo_array),max(dos.VovVo_array),n) linspace(min(dos.VovVo_array),max(dos.VovVo_array),n)];
Z_data = [dos.Evib_array  zeros(1,n) zeros(1,n)];
% figure(1)
% plot3(X_data,Y_data,Z_data,'ko');
% xlabel('X values')
% ylabel('Y values')
% zlabel('z values')
% grid on

% spline evaluation range (Should always be greater than data range)
X_range = [min(X_data) max(X_data)];
Y_range = [min(Y_data) max(Y_data)]; 

dZ = 1;

% chose number of control points for P and T and set the grid of control points

nxc=6;  % adjust these to chnage the numbers of control points
nyc=7;  % adjust these to chnage the numbers of control points
Xc=linspace(-1,X_range(2),nxc);
Xc=[0 logspace(-1,2,50) 150 200 300 700 1000 2000];
Yc=linspace(Y_range(1),Y_range(2),nyc);
options.Xc={Xc,Yc};

% The values given below were manually adjusted to achieve a reduced chi square of about 1
% change these numbers to see the effect on misfits and smoothness of the resulting surface
% smaller values for less misfit and less "smoothness", larger values for greater misfit and smoothness
options.mdrv = [2 2]; % derivative for which regularization is applied on X and Y (default 2nd)
options.lam=[90e2 2e1];  % intensity of regularization
options.ordr=[5 5]; % order of the spline


% RegFac is a vector of integers for the regularization oversampling (number of regularization pts between each control point)
%           if RegFac is not included in options, the default  is 1 for all dimensions     
%           if RegFac is included but empty, do a least square fit without regularization
%                (likely to fail unless data coverage is adequate)  
options.RegFac = [5 1];     

% all other options are allowed to be the default values

% since data are scattered, input X is a matrix with as many columns as
% independent variables
sp_E = spdft([X_data(:) Y_data(:)],Z_data(:),dZ(:),options);

%plot E
figure (1)
subplot(3,2,1)
fnplt(sp_E)
hold on
plot3(X_data,Y_data,Z_data,'ko','MarkerFaceColor','k')
hold off
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('Evib (J/kg)')

%plot residuals for the spline fit to the thermal energy 
subplot(3,2,2)
plot3(X_data,Y_data,(Z_data-fnval(sp_E,[X_data;Y_data]))./Z_data,'o','MarkerFaceColor','k')
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('Fractional residual')

%% Construct spline for isochoric heat capacity, entropy, plus ion-thermal 
   %contributions to the pressure, Helmholtz energy, and Gibbs energy

%Cv = temperature derivative of internal energy spline (J/kg/K)
sp_Cv = fnder(sp_E,[1 0]);
subplot(3,2,3)
fnplt(sp_Cv)
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('Cv (J/kg/K)')

%calculate Cv divided by T for every T,V point and fit a spline 
n = size(dos.T_array,2);
for i = 1:n
    T = dos.T_array(i);
    VovVo = dos.VovVo_array(i);
    if T == 0
        Cv_over_T(i) = 0;
    else
        Cv_over_T(i) = fnval(sp_Cv,[T;VovVo])/T;
    end
end
%fit a spline to the Cv/T points
X_data = [dos.T_array];
Y_data = [dos.VovVo_array];
Z_data = [Cv_over_T];
sp_Cv_over_T = spdft([X_data(:) Y_data(:)],Z_data(:),dZ(:),options);

%using the spline for Cv/T, determine the entropy S (J/kg/K) by
%integrating Cv/T with respect to temperature
sp_Cv_over_T_intT = fnder(sp_Cv_over_T,[-1 0]);
for i = 1:n
    T = dos.T_array(i);
    VovVo = dos.VovVo_array(i);
    dos.S_array(i) = fnval(sp_Cv_over_T_intT,[T;VovVo]);
end
%fit a spline to these discrete entropy points
Z_data = [dos.S_array];
sp_S = spdft([X_data(:) Y_data(:)],Z_data(:),dZ(:),options);
%plot
subplot(3,2,4)
fnplt(sp_S)
hold on
plot3(X_data,Y_data,Z_data,'ko','MarkerFaceColor','k')
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('S (J/kg/K)')

%using the spline for S, determine the thermal pressure Ptherm (GPa)
%sp_dSdVovVo_intT = spline representing temperature integral of the strain derivative 
%of entropy (since strain = VovVo is unitless, this integral has units of J/kg)
sp_dSdVovVo_intT = fnder(sp_S,[-1 1]);
%calculate Ptherm for every T,V point (make sure to convert units to GPa)
for i = 1:n
    T = dos.T_array(i);
    VovVo = dos.VovVo_array(i);
    dos.Ptherm_array(i) = fnval(sp_dSdVovVo_intT,[T;VovVo])/Vo*mw*1.e-6;
end
%fit a spline to these discrete points
Z_data = [dos.Ptherm_array];
sp_Ptherm = spdft([X_data(:) Y_data(:)],Z_data(:),dZ(:),options);
fnval(sp_Ptherm, [300; 0.7])
%plot
subplot(3,2,5)
fnplt(sp_Ptherm)
hold on
plot3(X_data,Y_data,Z_data,'ko','MarkerFaceColor','k')
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('P_{therm} (GPa)')

%compute the thermal contributions to the Helmholtz and Gibbs energies (J/kg)
for i = 1:n
    T = dos.T_array(i);
    VovVo = dos.VovVo_array(i);
    dos.Ftherm_array(i) = dos.Evib_array(i) - T*dos.S_array(i);
    dos.Gtherm_array(i) = dos.Ftherm_array(i) + ...
        dos.Ptherm_array(i)*Vo*1.e6/mw*VovVo;
end
%fit a spline to these discrete points
F_data = [dos.Ftherm_array];
G_data = [dos.Gtherm_array];
sp_Ftherm = spdft([X_data(:) Y_data(:)],F_data(:),dZ(:),options);
sp_Gtherm = spdft([X_data(:) Y_data(:)],G_data(:),dZ(:),options);
%plot
subplot(3,2,6)
fnplt(sp_Ftherm)
hold on
plot3(X_data,Y_data,F_data,'o','MarkerFaceColor','b')
fnplt(sp_Gtherm)
plot3(X_data,Y_data,G_data,'o','MarkerFaceColor','r')
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('Free energies (J/kg)')
