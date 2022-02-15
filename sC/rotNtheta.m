function out = rotNtheta(n, theta)
% ROTNTHETA(n, theta) returns the rotation matrix OUT associated to the 
%   general rotation of angle THETA around axis N.

    % TODO: why element-wise mult?
    out = n * n.' + (eye(3) - n * n.') .* cos(theta) + hat(n) .* sin(theta);

end