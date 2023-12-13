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
% total length). This allows for more simplified calculations. The
% following section defines the length, width, and frontal area of my kite.

W = 20; %length of beam (in)
L = 20;

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
V_air = 5.6; %m/s
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
[~,~,r1v,r2v,BiV] = bridlept(L,CoG,CoP,Bi,alpha); %calc reaction vectors

% plot bridal point, r1, r2 to check
bridleplotting(L,CoG,CoP,r1v,r2v,BiV);

%check for moment eq
[moment,~,~,~,~,~] = moments(A,alpha,V_air,rho_air,r1v,r2v,m);
%% Solving for Moments: The Moments (you have been waiting for)
BiT = linspace(.05,.06,100);
r1vT = zeros(100,3);
r2vT = zeros(100,3);
momT = zeros(100,3);

for i  = 1:length(BiT)
    [~,~,r1vT(i,:),r2vT(i,:)] = bridlept(L,CoG,CoP,BiT(i),alpha);
    [momT(i,:),~,~,~] = moments(A,alpha,V_air,rho_air,r1vT(i,:),r2vT(i,:),m);
end

momT = momT(:,3);
[~,idx] = min(abs(momT));
figure();plot(BiT,momT); hold on
plot(BiT(idx),momT(idx),'.','color','red','MarkerSize',12)
xlabel('Sampled Values of Bridle Point Offset')
ylabel('Torque (N-m)')
hold off

Bi = BiT(idx);
%% Final Kite Design
% Optimized for 5.5 m/s wind speed and a 12° angle of attack

[~,~,r1v,r2v,BiV] = bridlept(L,CoG,CoP,Bi,alpha);
[torque,F_d,F_l,F_g,F_bx,F_by] = moments(A,alpha,V_air,rho_air,r1v,r2v,m);
bridleplotting(L,CoG,CoP,r1v,r2v,BiV)

%vectorize forces to visualize the free body diagram

F_G_v = [0; F_g];
F_D_v = [F_d; 0];
F_L_v = [0; F_l];
F_B_v = [F_bx; F_by];
%% FBD
%fbd plot; looks bad, won't be trying to use matlab for this again.
R = [cos(-alpha), -sin(-alpha);
     sin(-alpha), cos(-alpha)];
F_G_v = R * F_G_v;
F_D_v = R * F_D_v;
F_L_v = R * F_L_v;
F_B_v = R * F_B_v;

figure();
    plot([0,L],[0,0],'Color','Black') %Kite side view
    %axis equal;
    hold on 
    plot(CoP,0,'.','Color','Red','MarkerSize',20) %CoP
    plot(CoG,0,'.','Color','Green','MarkerSize',20) %CoG
    plot(BiV(1),BiV(2),'.','Color','Blue','MarkerSize',20) %bridle pt
    %plot([BiX,CoP],[BiY,0],'--','Color','Red') %r1
    %plot([BiX,CoG],[BiY,0],'--','Color','Green') %r2
    quiver(BiV(1),BiV(2),r1v(1),r1v(2),'--','Color','Red')
    quiver(BiV(1),BiV(2),r2v(1),r2v(2),'--','Color','Green')
    quiver(CoG, 0, F_G_v(1),F_G_v(2)) %gravity force
    quiver(CoP, 0, F_D_v(1),F_L_v(2)) %pressure force
    %quiver(CoP, 0, F_L_v(1),F_L_v(2)) %lift force
    quiver(BiV(1),BiV(2),F_B_v(2),F_B_v(1)) %tension force
    legend('','CoP','CoG','Bridle Point','','','F_{grav}','F_{pressure}','F_{tension}','Location','Northwest')
    xlabel('Length (m)'); ylabel('Height (m)')
hold off

%% Plotting L_b
figure()
plot([0,L],[0,0],'color','blue'); hold on; axis equal
plot([0,L],[-L*sin(alpha),0],'--','color','black')
plot(CoP, 0, '.','Color','Green','Markersize', 20)
plot(CoG, 0, '.','Color','Red','Markersize', 20)
legend('','L_b','CoP','CoG')
xlabel('Length (m)')
ylabel('Height (m)')
hold off
%% dont forget about fsolve()
% I did! Oop