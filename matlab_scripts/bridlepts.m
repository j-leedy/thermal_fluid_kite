function [r1v,r2v] = bridlepts(alpha,x,L,d)
%BRIDAL POINTS calculates the bridal point and reaction vectors given an
%offset
%   This function is a helper function for a parameter sweep of the bridal
%   point, and therefore the radii r1,r2, that govern the moments acting on
%   a traditional diamond-style kite.
    %   INPUTS:
    %       alpha = angle of attack in radians
    %       x     = offset of bridal point with respect to end of kite
    %       L     = length of kite
    %       d     = distance between CoP and CoG on kite
    %
    %   OUTPUTS:
    %       r1v   = 3-D vector form of reaction vector btwn CoP and bridle
    %               point
    %       r2v   = 3-D vector btwn CoG and bridle point


x = L - x; % defining bridal coordinate along length of kite

% magnitude of r1
a = ((3*L)/4) - x; %length from bridal pt to CoP
r1 = a * tan(alpha); %magnitude of r1

% magnitude of r2
omega = atan(d/r1); %angle between r1 and r2
r2 = d/(sin(omega)); %magnitude of r2

% vectorize r1,r2
theta = (pi/2) - alpha; %angle between r1 and bridle cord
phi = theta - omega; %angle between r2 and bridle cord

r1v = [r1*cos(theta); 
       r1*sin(theta);
       0]; %3D vector form of r1
r2v = [r2*cos(phi);
       r2*sin(phi);
       0];


end