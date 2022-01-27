function H = GlobalConstraintMatrix(list)
% GLOBALCONSTRAINTMATRIX(list) returns the global constraint matrix H for 
%   the supplied list.
%   The structure of each element in LIST must follow the structure of 
%   ConstraintMatrix(). TODO: specify here too?

    consmats = [];
    for i = [1, size(list,1)]
        consmats = [consmats, ConstraintMatrix(list(i))];
    end

    H = diag(consmats);
end