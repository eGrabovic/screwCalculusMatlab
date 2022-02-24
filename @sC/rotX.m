function R = rotX(alpha)
% ROTX(alpha) computes the rotation matrix R related to a rotation of angle 
%   ALPHA around the x-axis.
%   Allows for symbolic compiutations.
            
    R = [[1,    0,             0];
         [0,    cos(alpha),    -sin(alpha)];
         [0,    sin(alpha),    cos(alpha)]];
    
end