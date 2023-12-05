function [eq_torque] = moments(A,alpha,V,rho,r1,r2,m)
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
F_drag = (C_d*rho*(V^2)*A_f)/2; %magnitude
F_lift = (C_l*rho*(V^2)*A)/2;
F_grav = -9.81*m;

%force vectors
F_g_v = [r2(1); r2(2)+F_grav; 0];
F_l_v = [r1(1); r1(2)+F_lift; 0];
F_d_b = [r1(1) + F_drag; r1(2); 0];

%sum of moments
eq_torque = cross(F_g_v,r2) + cross(F_l_v,r1) + cross(F_d_b,r1);

end

