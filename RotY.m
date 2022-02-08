function mat = RotY(theta, deg_bool)
% ROTY(theta) describes a rotation of the Euler angle THETA around
%   the Y-axis. If DEG_BOOL is true, the input angle will be considereed as
%   in degree.
%   Note that the output matrix is always with angles in radiants.
%
%   Examples
%       RotY(pi) or RotY(180, true)
%   --> [-1   0    0
%         0   1    0
%         0   0   -1]
%
%   Input
%       theta:      angle of rotation [rad or deg]
%       deg_bool:   boolean to consider input as degree
%   Output
%       mat:        3x3 rotation matrix
%

    if exist('deg_bool', 'var') && (deg_bool == true)
        th = deg2rad(theta);
    else % deg_bool does not exist or it's not a boolean 'true'
        th = theta;
    end

    mat = [cos(th),     0,  sin(th);...
           0,           1,  0;...
           -sin(th),    0,  cos(th)];

end