function [r1,r2,r1v,r2v] = bridlept(L,CoG,CoP,Bi,alpha)
%BRIDLEPT calculates a kite's reaction vectors
%   This function will calculate the reaction vectors pointing to the CoP
%   and CoG from a experimental value of bridal point. It will also output
%   the magnitude of those vectors.
%   INPUTS:
%       L     = Length of kite (m)
%       CoP   = X-component of center of gravity 
%       CoG   = X-component of center of mass
%       Bi    = X-component of proposed bridle point with respect to kite origin
%       alpha = Angle of attack (rad)
%   OUTPUTS:
%       r1    = Magnitude of r1 (reaction vector to CoP)
%       r2    = Magnitude of r2 (reaction vector to CoG)
%       r1v   = 3D vector representation of r1 with Bi as origin
%       r2V   = 3D vector representation of r2 with Bi as origin

%L_b = L*cos(alpha); %projection of kite L along angle of attack
a1 = L - CoG; %distance to CoG from near end of kite
a2 = a1 - CoP; %distance to CoP from same

% Calculate y-component of reaction vectors
Y_g = a1*sin(alpha); 
Y_p = a2*sin(alpha); 

% Calculate x-component of reaction vectors
X_g = a1*cos(alpha) - Bi; 
X_p = a2*cos(alpha) - Bi; 

% Creating Reaction Vectors
r1v = [X_p,Y_p,0];
r2v = [X_g,Y_g,0];

% Calculate Magnitude
theta = atan(Y_p/X_p); %angle between r1 and projected length
phi = atan(Y_g/X_g); %angle between r2 and same

r1 = sin(theta)/Y_p;
r2 = sin(phi)/Y_g;
end