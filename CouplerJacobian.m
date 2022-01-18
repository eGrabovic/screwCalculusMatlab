function res = CouplerJacobian(fstring, pars, constraintlist, pointlist)
    fFKin  = ToExpression(fstring + "ToMat");
    fJac   = ToExpression(fstring + "ToSpatialJac");
    
    H  = GlobalConstraintMatrix(constraintlist);
    G = [];
    
    for i = [1, size(pointlist,1)]
        G  = [G, GlobalGraspMatrix(fFKin(pars) * (pointlist(i)))];
    end
    GT = Transpose(G);
    
    Jo = ObjectJac(fJac, pars); 

    res = H * GT * Jo;
end