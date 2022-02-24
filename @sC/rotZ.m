function R = rotZ(alpha)
% ROTZ(alpha) computes the rotation matrix R related to a rotation of angle 
%   ALPHA around the z-axis.
%   Allows for symbolic compiutations.

    R = [[cos(alpha),   -sin(alpha),    0];
         [sin(alpha),   cos(alpha),     0];
         [0,            0,              1]];
    
end