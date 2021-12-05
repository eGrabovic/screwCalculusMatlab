function mat = RotZ(theta, deg_bool)
% ROTZ(theta) describes a rotation of the Euler angle THETA around
%   the Z-axis. If DEG_BOOL is true, the input angle will be considereed as
%   in degree.
%   Note that the output matrix is always with angles in radiants.
%
%   Examples
%       RotZ(pi) or RotZ(180, true)
%   --> [-1     0   0
%         0    -1   0
%         0     0   1]
%
%   Input
%       theta:      angle of rotation [rad or deg]
%       deg_bool:   boolean to consider input as degree
%   Output
%       mat:        3x3 rotation matrix
%


    if(nargin == 2) % 'degree' specified
        assert(class(deg_bool) == "logical", ...
           "Only booleans are allowed as second argument");
        if deg_bool
            th = deg2rad(theta);
        else
            th = theta;
        end
    else
        th = theta;
    end

    mat = [cos(th), -sin(th),   0;...
           sin(th), cos(th),    0;...
           0,       0,          1];
       
end