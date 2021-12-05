function rpy_angles = MatToRPY(rotmat, deg_bool)
% MATTORPY(rotmat, deg_bool) extracts the RPY_ANGLES from the ZYX Rotation 
%   matrix ROTMAT. If DEG_BOOL is true, the angles are considered to be in 
%   degrees.
%
%   Examples...
%
%   Input
%       rotmat:     3x3 rotation matrix
%       deg_bool:   boolean to indicate that the angles are in degrees
%   Output
%       rpy_angles: 1x3 RPY angles [rad]
%
    if(nargin == 1) % nothing specified
        is_deg = false;
    else
        assert(class(deg_bool) == "logical", ...
           "Only booleans are allowed as second argument");
        is_deg = deg_bool;
    end

    rpy_angles = MatToEulZYX(rotmat, is_deg);
end