%% Kite Project 
% Thermal Fluid Systems FA23
% 
% By Joe Leedy
%%
close all
clear 
clc

%% Kite Shape Parameters
% Kite is a simple diamond design with a cross beam at the COP (1/4 the
% total length). This allows for more simplified complications. The
% following section defines the length, width, and frontal area of my kite.

W = 29; %length of beam (in)
L = 29;

W = W * 0.0254; %convert to meters
L = L * 0.0254;

A = (W*L)/2; % m^2 (surface area)

% kite mass
rods = 26 / 1000; %mass in kg of dowel cross bars

rho_paper = 1.15e3; %kg/m3
thck_paper = .1 / 1000; %m
paper = rho_paper*thck_paper*A; %mass of paper in kg
m = paper + rods; %mass of kite assuming tape + string are negligible

% angle of attack
alpha = deg2rad(12);
%% Air Properties

t_air = ((40-32)/1.8); % °C
hum_air = 40; % relative humidity %

%credit to @sjfitz on github for this funciton
[rho_air,mu_air] = AirProperties(t_air,[],hum_air); 

mu_air = mu_air * 0.1019; %conversion factor to kg/m-s

t_air = t_air + 273; %convert to K
V_air = 4; %m/s
%% Solving for Center of Pressure and Center of Gravity
kiteshape = [0 W/2 0;0 L/4 L];
% code for generating a peicewise function of the kite
m1 = (W/2)/(kiteshape(2,2)-kiteshape(2,1)); %eq. 1
b1 = 0;
m2 = -(W/2)/(kiteshape(2,3)-kiteshape(2,2)); %eq. 2
b2 = -m2*L;

%code for calculating center of gravity
Ad = @(x) x.*m1 + b1; %functions for denominator of CoG formula
Bd = @(x) x.*m2 + b2;

An = @(x) x.*(x.*m1 + b1); %functions for numerator 
Bn = @(x) x.*(x.*m2 + b2);

CoG = (integral(An, kiteshape(2,1),kiteshape(2,2))+ ...
       integral(Bn,kiteshape(2,2),kiteshape(2,3)))/ ...
       (integral(Ad, kiteshape(2,1),kiteshape(2,2))+ ...
       integral(Bd,kiteshape(2,2),kiteshape(2,3)));

CoP = L/4; %CoP is one fourth the relevant chord length 

plotkite(kiteshape,CoP,CoG);
%% Solving for Moments: Bridal Point Calculation

Bi = .1; %set initial bridle point
[r1,r2,r1v,r2v] = bridlept(L,CoG,CoP,Bi,alpha);

% plot bridal point, r1, r2 to check
bridleplotting(L,CoG,CoP,r1v);

%check for moment eq
moments(A,alpha,V_air,rho_air,r1v,r2v,m)
%% Solving for Moments: The Moments (you have been waiting for)
%{
% x and y components of r1, r2
r1_v = [r1*cos(theta);
        r1*sin(theta);
        0]; %r1 in vector form using bridle point as origin
r2_v = [r2*cos(phi);
        r2*sin(phi);
        0]; %r2 in vector form

% parameter sweep on bridal point

x_trial = linspace(1,13,100) .* 0.0254; %experimental values for bridal point
x2_tmp = 6 * 0.0254; 
r1vt = zeros(3,1);
r2vt = zeros(3,1);
moment_t = zeros(3,length(x_trial),1);

%this doesn't work! evaluate with emily in class tmmr
for i = 1:length(x_trial)
    [r1vt,r2vt] = bridlepts(alpha,x_trial(i),L,d);
    moment_t(:,i) = moments(A,alpha,V_air,rho_air,r1vt,r2vt,m);
end
torques = moment_t(3,:);
[~,ideal_idx] = sort(torques,'ascend');
%}
%bridleplot(L,x_trial(100),x2_tmp,CoP,CoG,norm(r1vt))
%% dont forget about fsolve()