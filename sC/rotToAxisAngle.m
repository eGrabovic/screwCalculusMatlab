function [n, theta] = rotToAxisAngle(R)
% ROTTOAXISANGLE(R) uses the rotation matrix R to return the pair 
%   [N, THETA] where:
%   - N is the axis around which the rotation takes place
%   - THETA is the angle of rotation.

    theta = acos((trace(R) - 1) / 2) ;

    n = 1/(2*sin(theta)) .* ...
        [R(3, 2) - R(2, 3);
         R(1, 3) - R(3, 1);
         R(2, 1) - R(1, 2)];

end