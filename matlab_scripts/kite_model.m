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

W = 8.5; %length of beam (in)
L = 11.5;

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

CoP = L/4;
plot(CoP,0, 'o',"Color",'red','MarkerSize',10)
legend('','', 'Center of Pressure')

hold off