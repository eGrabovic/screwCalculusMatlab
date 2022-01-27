function fhJac = FreeHandJacobian(list)
% FREEHANDJACOBIAN(list) returns the global free hand Jacobian. This is the
%   block-diag Jacobians for all the fingers, where each one is described 
%   through its DH table, its joint variables and its initial pose w.r.t. 
%   to the base frame. All this is contained in a single element of the
%   supplied LIST.

    % each element of list must be [DHtable_i, vars_i, gp0i]
    Jacmats = [];
    for i = [1, size(list,1)]
        curr_list = list(i);
        fun = list(1);
        Jacmats = [Jacmats, DHJacobBase(fun(curr_list(2)), curr_list(3))];
    end

    fhJac   = diag(Jacmats);
end