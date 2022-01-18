function fhJack = FreeHandJacobian(list)
    % each element of list must be [DHtable_i, vars_i, gp0i]
    Jacmats = [];
    for i = [1, size(list,1)]
        curr_list = list(i);
        fun = list(1);
        Jacmats = [Jacmats, DHJacobBase(fun(curr_list(2)), curr_list(3))];
    end

    fhJac   = diag(Jacmats);
end