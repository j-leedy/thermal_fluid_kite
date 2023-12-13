function [r1,r2,r1v,r2v,BiV] = bridlept(L,CoG,CoP,Bi,alpha)
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

Lb = L/cos(alpha); %Projection of kite along AoA
c = Lb-Bi; %Distance between far end of kite and bridle point
b = c*cos(alpha); %x-distance between bridle point and kite origin
a = c*sin(alpha); %y-distance between and of kite and projected bridle pt

% Calculate y-component of reaction vectors
Y_p = a; %y-component same for both
Y_g = a;

% Calculate x-component of reaction vectors
BiX = L- b; %x component of projected bridle pt w ref to kite origin
X_p = CoP - BiX; 
X_g = CoG - BiX;

% Vector form of bridle pt (wrt kite origin)
BiV = [BiX; -a; 0];

% Creating Reaction Vectors
r1v = [X_p,Y_p,0];
r2v = [X_g,Y_g,0];

% Calculate Magnitude
r1 = norm(r1v);
r2 = norm(r2v);
end