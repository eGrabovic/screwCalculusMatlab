function res = PfaffianMatrix(fingerlist, fstring, pars, constraintlist, pointlist)
    A11 = HandJacobian(fingerlist, constraintlist);
    
    A12 = - CouplerJacobian(fstring, pars, constraintlist, pointlist);
    
    res = [A11; A12];
end