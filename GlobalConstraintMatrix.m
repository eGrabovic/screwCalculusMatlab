function H = GlobalConstraintMatrix(list)
    consmats = [];
    for i = [1, size(list,1)]
        consmats = [consmats, ConstraintMatrix(list(i))];
    end

    H = diag(consmats);
end