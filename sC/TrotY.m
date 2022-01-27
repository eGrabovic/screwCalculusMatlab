function out = TrotY(x)
% TODO


if isa(x, 'ADvar')
    T = rotY(x.val);
    Der = hat([0;1;0])*T.*x.der;
    Der = [Der, [0;0;0];0 0 0 0];
    out = ADvar(TrotY(x.val), Der);
    return
end


out = [cos(x), 0, sin(x),0;...
       0,      1,     0, 0;...
      -sin(x), 0, cos(x),0;...
       0,      0,      0,1];

end