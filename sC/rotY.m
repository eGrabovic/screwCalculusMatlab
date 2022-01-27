function R = rotY(alpha)
% ROTY(alpha) computes the rotation matrix R along the Y-axis.

if isa(alpha, 'ADvar')
    T = rotY(alpha.val);
    Der = hat([0;1;0])*T.*alpha.der;
    R = ADvar(T, Der);
    return
end
R = [cos(alpha) 0 sin(alpha);0 1 0; -sin(alpha) 0 cos(alpha)];
end