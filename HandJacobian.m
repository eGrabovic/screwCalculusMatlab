function hJac = HandJacobian(fingerlist, constraintlist)
% HANDJACOBIAN(fingerlist, constraintlist) returns the H*J matrix HJAC 
%   given the finger list FINGELIST as in FreeHandJacobian() and CONSTRAINTLIST
%   as in GlobalConstraintMatrix().
%   TODO: maybe copy the info here too?

    fhJac = FreeHandJacobian(fingerlist);
    H = GlobalConstraintMatrix(constraintlist);
    
    hJac  = H * fhJac;
end