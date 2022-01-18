function hJac = HandJacobian(fingerlist, constraintlist)
    fhJac = FreeHandJacobian(fingerlist);
    H = GlobalConstraintMatrix(constraintlist);
    
    hJac  = H * fhJac;
end