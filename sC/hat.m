function X = hat(vec)
% HAT(vec) transforms a R3 column vector VEC in its 'hat' R3x3
%   antisymmetric matrix form X, or transforms a R3 column twist vector VEC
%   in a homogeneous 'hat' R4x4 matrix X.
%   Allows for ADvar class argument.

    % Manage ADvar argument
    if isa(vec, 'ADvar')
        X = ADvar(hat(vec.val), hat(vec.der));
        return
    end
    
    % R3x3 column vector case
    if size(vec, 1) == 3
        X = [[0,         -vec(3),      vec(2)]
             [vec(3),     0,           -vec(1)]
             [-vec(2),    vec(1),      0]];
        return
    end
    
    % TODO: remove from here
    % R3 (6-element) twist vector case
    X = [[hat([vec(4); vec(5); vec(6)]),  vec(1:3)]
         [0,   0,   0,                    0]];

end