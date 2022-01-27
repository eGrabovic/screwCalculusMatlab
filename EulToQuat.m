function quat = EulToQuat(EU_angles, type, deg_bool)
% EULTOQUAT(EU_angles, type, deg_bool) returns the unit quaternion QUAT
%   which describes the rotation related to the specified angles EU_ANGLES 
%   and parametrization TYPE.
%   Valid parametrizations consist of a combination of 'X', 'Y', 'Z'
%   (uppercase).
%   If a boolean DEG_BOOL is specified as 'true', input angles will be
%   considered as in degrees.
%   Note that the output matrix is always with angles in radiants.

    quat = MatToQuat(EulToMat(EU_angles,type, deg_bool));
end