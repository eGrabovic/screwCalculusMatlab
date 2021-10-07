function mat = RotZ(theta)
% ROTZ(theta) describes a rotation of the Euler angle THETA around
%   the Z-axis.
%
%   mat = ROTZ(theta) return the rotation matrix.
%
%   theta - angle of rotation [rad]
%

    mat = [cos(theta),-sin(theta),0;...
           sin(theta),cos(theta),0;...
           0,0,1];
       
end