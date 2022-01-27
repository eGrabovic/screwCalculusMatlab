function R = rotZ(alpha)
% ROTZ(alpha) computes the rotation matrix R along the Z-axis.

if isa(alpha, 'ADvar')
    T = rotZ(alpha.val);
    Der = hat([0;0;1])*T.*alpha.der;
    R = ADvar(T, Der);
    return
end

R = [cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1];
end