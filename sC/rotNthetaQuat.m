function out = rotNthetaQuat(n, theta)
% ROTNTHETAQUAT returns the rotation OUT associated to the general rotation
%   of angle THETA around the axis N. Internally, it uses a quaternion
%   representation.

    p0 = cos(theta./2);
    p1 = sin(theta./2) .* n(1);
    p2 = sin(theta./2) .* n(2);
    p3 = sin(theta./2) .* n(3);

    out = NaN(3); % 3x3 NaN matrix

    % Build rotation matrix
    out(1,1) = p0.^2 + p1.^2 - 0.5;
    out(2,1) = p1.*p2 - p0.*p3;
    out(3,1) = p1.*p3 + p0.*p2;

    out(1,2) = p1.*p2 + p0.*p3;
    out(2,2) = p0.^2 + p2.^2 - 0.5;
    out(3,2) = p2.*p3 - p0.*p1;

    out(1,3) = p1.*p3 - p0.*p2;
    out(2,3) = p2.*p3 + p0.*p1;
    out(3,3) = p0.^2 + p3.^2 - 0.5;

    out = 2 .* out.';
end