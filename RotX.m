function mat = RotX(theta, deg_bool)
% ROTX(theta) describes a rotation of the Euler angle THETA around
%   the X-axis. If DEG_BOOL is true, the input angle will be considereed as
%   in degree.
%   Note that the output matrix is always with angles in radiants.
%
%   Examples
%       RotX(pi) or RotX(180, true)
%   --> [1   0   0
%        0  -1   0
%        0   0  -1]
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

    mat = [1,   0,          0;...
           0,   cos(th),    -sin(th);...
           0,   sin(th),    cos(th)];

end