function X = hat(screwObj, v)
% HAT(screwObj, v) transforms a R3 column vector V in its 'hat' R3x3
%   antisymmetric matrix form X, or transforms a R3 column twist vector V 
%   in a homogeneous 'hat' R4x4 matrix X.
%   Allows for ADvar class argument.

    % Manage ADvar argument
    if isa(vec, 'ADvar')
        X = ADvar(hat(vec.val), hat(vec.der));
        return
    end

    if all(size(v) ~= [3,1]) || all(size(v) ~= [6,1])
        error('v must be a 3x1 matrix or a 6-element vector')
    end
    
    % R3x3 column vector case
    if all(size(v) == [3,1])
        X = [[ 0,        -v(3),    v(2)];
             [ v(3),      0,      -v(1)];
             [-v(2),      v(1),    0]];
    end
    
    % TODO: remove from here?
    % R3 (6-element) twist vector case
    if all(size(v) == [6,1])
        vel = [v(1); v(2); v(3)];
        omeg = [v(4); v(5); v (6)];

        X = [[screwObj.hat(omeg),   vel];
             [zeros(1, 3),          0]];
    end

end