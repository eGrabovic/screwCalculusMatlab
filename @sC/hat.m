function X = hat(v)
% HAT(vec) transforms a R3 column vector V in its 'hat' R3x3
%   antisymmetric matrix form X, or transforms a R3 column twist vector V 
%   in a homogeneous 'hat' R4x4 matrix X.  

    if all(size(v) ~= [3,1]) || all(size(v) ~= [6,1])
        error('v must be a 3x1 matrix or a 6-element vector')
    end
    
    if all(size(v) == [3,1])
        X = [[0,        -v(3),  v(2)];
             [v(3),     0,      -v(1)];
             [-v(2),    v(1),   0]];
    end
    
    if all(size(v) == [6,1])
        vel = [v(1); v(2); v(3)];
        omeg = [v(4); v(5);v (6)];
        X = [[sC.hat(omeg), vel];
             [0, 0, 0,      0]];
    end

end