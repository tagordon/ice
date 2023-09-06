%% Import data
 close all
 clear all

%addpath('./Thermo_LBF_Tools/');
tic
%high_pressure_ice_DOS;
toc

addpath('../')
%iceVII_DOS_family

%save dos.mat dos
%%
close all
clear all
load dos.mat
tic
n = 50;
% X_data = [dos.T_array zeros(1,n) ones(1,n)*0.02 ones(1,n)*2 ones(1,n)*3  ones(1,n)*4 ones(1,n)*10.0];
% Y_data = [dos.VovVo  linspace(0.1,1,n) linspace(0.1,1,n) linspace(0.1,1,n)  linspace(0.1,1,n) linspace(0.1,1,n) linspace(0.2,1,n)];
% Z_data = [dos.Evib_array  zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n)];
 X_data = [dos.T_array zeros(1,n) ones(1,n) ones(1,n)*2];
 Y_data = [dos.VovVo linspace(0.15,1,n) linspace(0.15,1,n) linspace(0.15,1,n)];
 Z_data = [dos.Evib_array   zeros(1,n) zeros(1,n) zeros(1,n)];
% Check input data

% X_data = [X_data linspace(0,300,40) ones(1,50) ];
% Y_data = [Y_data linspace(0.35,0.1,40) linspace(1,0.7,50) ];
% Z_data = [Z_data zeros(1,40) zeros(1,50)];


% id=find(dos.T_array==100);
% p = polyfit([dos.VovVo(id) 0 0 0 0],[dos.Evib_array(id) 0 0 0 0],3);
% V=linspace(0.1,0.35,n);
% E=p(4)+p(3).*V + p(2).*V.^2 + p(1).*V.^3;
% X_data = [X_data ones(1,n)*100];
% Y_data = [Y_data V ];
% Z_data = [Z_data  E ];




%id=find(dos.T_array==150);
%p = polyfit([dos.VovVo(id) 0 0 0 0],[dos.Evib_array(id) 0 0 0 0],3);
%V=linspace(0.1,0.35,n);
%E=p(4)+p(3).*V + p(2).*V.^2 + p(1).*V.^3;
%X_data = [X_data ones(1,n)*150];
%Y_data = [Y_data V ];
%Z_data = [Z_data  E ];

% id=find(dos.T_array==100);
% p = polyfit([dos.VovVo(id) 0 0 0 0],[dos.Evib_array(id) 0 0 0 0],3);
% V=linspace(0.1,0.35,n);
% E=p(4)+p(3).*V + p(2).*V.^2 + p(1).*V.^3;
% X_data = [X_data ones(1,n)*150];
% Y_data = [Y_data V ];
% Z_data = [Z_data  E ];

% id=find(dos.T_array==50);
% p = polyfit([dos.VovVo(id) 0 0 0 0],[dos.Evib_array(id) 0 0 0 0],3);
% V=linspace(0.1,0.35,n);
% E=p(4)+p(3).*V + p(2).*V.^2 + p(1).*V.^3;
% X_data = [X_data ones(1,n)*150];
% Y_data = [Y_data V ];
% Z_data = [Z_data  E ];
% 
% 
% 
% 
% id=find(dos.T_array==1300);
% p = polyfit([dos.VovVo(id)],[dos.Evib_array(id)],3);
% V=linspace(0.1,0.35,n);
% E=p(4)+p(3).*V + p(2).*V.^2 + p(1).*V.^3;
% X_data = [X_data ones(1,n)*1300];
% Y_data = [Y_data V ];
% Z_data = [Z_data  E ];

figure(2)
plot3(X_data,Y_data,Z_data,'ko');
xlabel('X values')
ylabel('Y values')
zlabel('z values')
grid on
toc

%% Evaluating spline with regularization
tic

% spline evaluation range (Should always be greater than data range)
X_range = [min(X_data) max(X_data)];
Y_range = [0.1 1]; 

dZ = 1;

% chose number of control points for P and T and set the grid of control points
%!!!! fewer knots
nxc=10;  % adjust these to chnage the numbers of control points 
nyc=10;  % adjust these to chnage the numbers of control points 
Xc=linspace(-1,X_range(2),nxc);
Xc=[-1 5 50 100 150 200 300 500 700 1000 1300 1500 1700 2000]; %!!!! play with knots
Xc = [-5, 0, 5, 20, 40, 60, 80, 100, 200, 300, 400, 500, 600, 800, 1000, 1400, 1800, 2000];
Yc=linspace(Y_range(1),Y_range(2),nyc);
Yc = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
options.Xc={Xc,Yc};



% The values given below were manually adjusted to achieve a reduced chi square of about 1
% change these numbers to see the effect on misfits and smoothness of the resulting surface
% smaller values for less misfit and less "smoothness", larger values for greater misfit and smoothness

options.mdrv = [4, 2]; % derivative for which regularization is applied on X and Y (default 2nd)

options.lam=[7e8 1];  % intensity of regularization

options.ordr=[7 6]; % order of the spline !!!! higher order


% RegFac is a vector of integers for the regularization oversampling (number of regularization pts between each control point)
%           if RegFac is not included in options, the default  is 1 for all dimensions     
%           if RegFac is included but empty, do a least square fit without regularization
%                (likely to fail unless data coverage is adequate)  
options.RegFac = [1 1];

% all other options are allowed to be the default values

% since data are scattered, input X is a matrix with as many columns as
% independent variables

sp=spdft([X_data(:) Y_data(:)],Z_data(:),dZ(:),options);


figure (1)
subplot(2,2,1)
fnplt(sp)
hold on
plot3(X_data,Y_data,Z_data,'ko','MarkerFaceColor','k')
hold off
xlabel('Temperature (K)')
ylabel('V/Vo')
zlabel('Evib (J/kg)')

sp_Cv= fnder(sp,[1 0]);

subplot(2,2,2)
plot3(X_data,Y_data,(Z_data-fnval(sp,[X_data ;Y_data])),'ko','MarkerFaceColor','k')
grid on
subplot(2,2,3)
fnplt(sp_Cv)

toc
%%
spEvib = sp;
save spEvib_VII.mat spEvib