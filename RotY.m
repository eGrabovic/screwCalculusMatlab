function mat = RotY(theta)
% ROTY(theta) describes a rotation of the Euler angle THETA around
%   the Y-axis.
%
%   mat = ROTY(theta) return the rotation matrix.
%
%   theta - angle of rotation [rad]
%

    mat = [cos(theta),  0,  sin(theta);...
           0,           1,  0;...
           -sin(theta), 0,  cos(theta)];

end