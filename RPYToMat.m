function mat = RPYToMat(angles, deg_bool)
% RPYTOMAT(angles, deg_bool) considers the ANGLES specified as Euler angles
%   used to build a ZYX Rotation matrix. If DEG_BOOL is true, the angles
%   are considered to be in degrees.
%
%   Examples
%       RPYToMat([0,0,pi/4]) or RPYToMat([0,0,45], true)
%   --> [1    0         0
%        0    0.7071   -0.7071
%        0    0.7071    0.7071]
%
%   Input
%       angles:     1x3 numerical array of angles
%       deg_bool:   boolean to indicate that the angles are in degrees
%   Output
%       mat:        3x3 Rotation matrix
%
    if(nargin == 1) % nothing specified
        is_deg = false;
    else
        assert(class(deg_bool) == "logical", ...
           "Only booleans are allowed as second argument");
        is_deg = deg_bool;
    end

    mat = EulToMat(angles,"ZYX", is_deg);
end