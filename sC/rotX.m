function R = rotX(alpha)

if isa(alpha, 'ADvar')
    T = rotX(alpha.val);
    Der = hat([1;0;0])*T.*alpha.der;
    R = ADvar(T, Der);
    return
end

R = [1 0 0;0 cos(alpha) -sin(alpha);0 sin(alpha) cos(alpha)];

end