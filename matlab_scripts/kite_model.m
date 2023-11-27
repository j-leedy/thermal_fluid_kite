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
%% Solving for CoP and CoG

kiteshape = [0 W/2 0;0 L/4 L]; %points to plot the kite shape
cline = [0 0;0 L];

%plot the kite and a dashed centerline
plot(kiteshape(2,:),kiteshape(1,:)); axis equal
hold on; plot(cline(2,:),cline(1,:), '--'); 

CoP = [L/4 ; 0]; %i'm not sure why this isn't currently working but its not!
plot(CoP, 'o',"Color",'black','MarkerSize',10)
hold off