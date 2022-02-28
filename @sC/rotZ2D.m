function R = rotZ2D(screwObj, alpha)
% ROTZ2D(alpha) computes the 2D rotation matrix R related to a rotation of 
%   angle ALPHA around the x-axis. 
%   Allows for symbolic computations.
    
    R = [[cos(alpha),   -sin(alpha)];
         [sin(alpha),   cos(alpha)]];
end