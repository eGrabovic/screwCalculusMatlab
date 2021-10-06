function R = rotY(alpha)
%
% funzione matrice di rotazione attorno asse Y per acconsentire
% il calcolo simbolico
%
if isa(alpha, 'ADvar')
    T = rotY(alpha.val);
    Der = hat([0;1;0])*T.*alpha.der;
    R = ADvar(T, Der);
    return
end
R = [cos(alpha) 0 sin(alpha);0 1 0; -sin(alpha) 0 cos(alpha)];
end