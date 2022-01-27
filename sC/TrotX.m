function out = TrotX(x)
% TODO

if isa(x, 'ADvar')
    T = rotX(x.val);
    Der = hat([1;0;0])*T.*x.der;
    Der = [Der, [0;0;0];0 0 0 0];
    out = ADvar(TrotX(x.val), Der);
    return
end

out = [1,      0,       0, 0;...
    0, cos(x), -sin(x),0;...
    0, sin(x),  cos(x),0;...
    0,0,0,1];

end