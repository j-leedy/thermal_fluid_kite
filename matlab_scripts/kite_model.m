% Intro Thermal-Fluid Systems
% Joe Leedy
% Kite Project 
%%
close all
clear all
clc

%% Kite Shape Parameters
% Kite is a simple diamond design with a cross beam
% at the COP (1/4 the total length)

W = 29; %length of beam (in)
L = 29;

W = W * 0.0254; %convert to meters
L = L * 0.0254;

A = (W*L)/2; % m^2 (surface area)
%% Air Properties

t_air = ((40-32)/1.8); % c
humidity = 40; % relative humidity %

[rho_air,mu_air] = AirProperties(t_air,[],humidity); %credit to @sjfitz on github for this funciton
mu_air = mu_air * 0.1019; %conversion factor to kg/m-s
%% Solving for CoP and CoG

kiteshape = [0 W/2 0;0 L/4 L]; %points to plot the kite shape
cline = [0 0;0 L];

%plot the kite and a dashed centerline
plot(kiteshape(2,:),kiteshape(1,:),'Color','blue'); axis equal
hold on; plot(cline(2,:),cline(1,:), '--', 'Color','black'); 

% code for generating a peicewise function of the kite
m1 = (W/2)/(kiteshape(2,2)-kiteshape(2,1)); %eq. 1
b1 = 0;
m2 = -(W/2)/(kiteshape(2,3)-kiteshape(2,2)); %eq. 2
b2 = -m2*L;

%code for calculating center of gravity
Ad = @(x) x.*m1 + b1;
Bd = @(x) x.*m2 + b2;

An = @(x) x.*(x.*m1 + b1);
Bn = @(x) x.*(x.*m2 + b2);

CoG = (integral(An, kiteshape(2,1),kiteshape(2,2))+ ...
       integral(Bn,kiteshape(2,2),kiteshape(2,3)))/ ...
       (integral(Ad, kiteshape(2,1),kiteshape(2,2))+ ...
       integral(Bd,kiteshape(2,2),kiteshape(2,3)));

%CoP is one fourth the relevant chord length 
CoP = L/4;
plot(CoP,0, 'o',"Color",'red','MarkerSize',10)
plot(CoG,0,'o','Color','green','MarkerSize',10)
legend('','', 'Center of Pressure','Center of Gravity')

hold off
%% dont forget about fsolve()
% for some reason matlab doesn't think the function exists?