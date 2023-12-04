function [outputArg1,outputArg2] = moments(A,alpha,V,rho,r1,r2,m)
% KITE_PARAMETER_TEST function to model performance of a simple kite
%   
%   Given inputs of surface area, angle of attack, wind velocity, ... , 
%       check for static equilibrium in x, y, and z; assume the system 
%       can be modelled as planar (2D).

%solving for coefficients
A_f = A*sin(alpha); %frontal area
C_d = 1.28 * sin(alpha); %coeff of drag
C_l = 2 * sin(alpha); %coeff of lift

%solving for forces
F_drag = (C_d*rho*(V^2)*A_f)/2;
F_lift = (C_l*rho*(V^2)*A)/2;
F_grav = -9.81*m;

end

