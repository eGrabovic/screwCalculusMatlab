function res = NullBasis(A, independent)
% NULLBASIS(A, independent) returns a matrix RES whose columns span the 
%   NullSpace of the input matrix A. The controls are the indexes of the 
%   components of q_dot in the list INDEPENDENT

% The whole expression A(q1,q2,...) must be supplied
    [m, n] = size(A);
    unit = eye(n - m);
    r = [1:n];   
    dependent = r(setdiff(1:end, independent));
    depsel = dependent(:, 1:n); %SelectionMatrixColumns(n, dependent);
    indepsel = independent(:,1:n); %SelectionMatrixColumns(n, independent);
    LHSmat =  A * depsel;
    RHSvec = -A * indepsel;
    res = inv(LHSmat) * RHSvec;
    
    for i = [1, size(independent, 2)]
        res = [res(:,1:independent(i)-1), unit(:,i), res(:, 1:independent(i)+1)];
    end
end