function or_err = MatToVect(R_act, R_des)
% MATTOVECT(R) returns the orientation error OR_ERR between the actual 
%   orientation related to the rotation matrix R_ACT and the desired one 
%   related to R_DES, such that 
%   or_err = r * sin(theta)
%   where r and theta are the components of the axis angle parametrization.
%   OR_ERR represents also the axis form of the skew-symm matrix 
%   R_SS = (1/2) * (R - R^T).
% TODO: how to add formulas?

    or_err = FramesToVect(R_act, R_des); 
end