function out = axisAngleToQuat(theta, n)
% AXISANGLETOQUAT returns the quaternion OUT associated with the rotation
%   of angle THETA around the general axis N.

    assert(isnumeric(theta), ...
           "Only numeric angle values are supported.");
    assert((max(size(n)) == 3) && (min(size(n)) == 1), ...
           "You must provide a 3-element array for the axis.")

    % Init
    out = zeros(4,1);

    % TODO: why element-wise product?
    out(1) = cos(theta/2);
    out(2) = sin(theta/2) .* n(1);
    out(3) = sin(theta/2) .* n(2);
    out(4) = sin(theta/2) .* n(3);
end