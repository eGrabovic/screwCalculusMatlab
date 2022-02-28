function R = rotY(screwObj, alpha)
% ROTY(alpha) computes the rotation matrix R related to a rotation of angle 
%   ALPHA around the y-axis.
%   Allows for symbolic compiutations.

    R = [[cos(alpha),   0,  sin(alpha)];
         [0,            1,  0];
         [-sin(alpha),  0,  cos(alpha)]];
end