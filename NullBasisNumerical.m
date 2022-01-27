function res = NullBasisNumerical(A, q0, independent)
% NULLBASISNUMERICAL(A, q0, independent) returns a matrix RES whose columns
%   span the NullSpace of the input matrix A. The controls are the indexes
%   of the components of q_dot in the list INDEPENDENT

% Only the head of the matrix has to be supplied
    
    if ismatrix(A(q0))
        Aq0 = A(q0);
    else
        Aq0 = A(:, q0);
    end
    
    [m, n] = size(Aq0);
    unit = eye(n - m);
    
    r = [1:n];
    dependent = r(setdiff(1:end, independent));
    depsel   = dependent(:,1:n); %SelectionMatrixColumns(n, dependent);
    indepsel = independent(:,1:n); %SelectionMatrixColumns(n, independent);
    LHSmat =  Aq0 * depsel;
    RHSvec = -Aq0 * indepsel;
    res = pinv(LHSmat) * RHSvec;
    
    for i = [1, size(independent,2)]
        res = [res(:,1:independent(i)-1), unit(:,i), res(:, 1:independent(i)+1)];
    end
end