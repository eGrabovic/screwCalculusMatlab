function X = hat(vec)
% HAT(vec) transforms a R3 column vector VEC in its 'hat' R3x3
% antisymmetric matrix form, or transforms a R3 column twist vector in a
% homogeneous 'hat' R4x4 matrix.         

        if isa(vec, 'ADvar')
            X = ADvar(hat(vec.val), hat(vec.der));
            return
        end
        
        if size(vec, 1) == 3
            X = [0 -vec(3) vec(2);...
                vec(3) 0 -vec(1);...
                -vec(2) vec(1) 0];
            return
        end
        
        X = [hat([vec(4);vec(5);vec(6)]), vec(1:3);0 0 0 0];


end