function R = rotZ(alpha)
%
% funzione matrice di rotazione attorno asse Z per acconsentire
% il calcolo simbolico
%

if isa(alpha, 'ADvar')
    T = rotZ(alpha.val);
    Der = hat([0;0;1])*T.*alpha.der;
    R = ADvar(T, Der);
    return
end

R = [cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1];
end