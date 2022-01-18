function [H, Adg] = ConstraintMatrix(constrainttype, axis, Rbe)
    H = eye(6);
    
    p = [0, 0, 0];
    RT = Rbe';
    
    g = RPToHomogeneous(RT, p);
    Adg = RigidAdjoint(g);
    
    if constrainttype == "S"
        H = H(1:3, 1:6); %Take(H, {1, 3}, {1, 6});
    end
    
    if constrainttype == "R"
        F = [[0; 0; 0], axis];
        H = NullSpace(F);
    end
    
    if constrainttype == "P"
        F = [axis, [0; 0; 0]];
        H = NullSpace(F);
    end
    
    if constrainttype == "C"
        H = eye(6);
    end
    
    if constrainttype == "PC"
        H = Take(H, {3});
    end
    
    if constrainttype == "PCWF"
        H = H(1:3,1:6); %Take(H, {1, 3}, {1, 6});
    end
    
    if constrainttype == "SF"
        H = H([1, 2, 3, 6], :);
    end
end