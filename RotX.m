function mat = RotX(theta)
% ROTX(theta) describes a rotation of the Euler angle THETA around
%   the X-axis.
%
%   mat = ROTX(theta) return the rotation matrix.
%
%   theta - angle of rotation [rad]
%

    mat = [1,0,0;...
           0,cos(theta),-sin(theta);...
           0,sin(theta),cos(theta)];

end