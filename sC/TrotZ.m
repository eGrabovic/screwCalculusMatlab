function out = TrotZ(alpha)
% TODO

if isa(alpha, 'ADvar')
    T = rotZ(alpha.val);
    Der = hat([0;0;1])*T.*alpha.der;
    Der = [Der, [0;0;0];0 0 0 0];
    out = ADvar(TrotZ(alpha.val), Der);
    return
end

out = [[cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1], [0;0;0];0,0,0,1];

end