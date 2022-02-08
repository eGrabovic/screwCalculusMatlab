function R = rotZ2D(alpha)
% ROTZ2D(alpha) returns the 2D rotation matrix R corresponding to a X-axis
%   rotation of angle ALPHA.
    
    R = [cos(alpha) -sin(alpha);
         sin(alpha)  cos(alpha)];
end