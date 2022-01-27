function mat = NDMatrix(S, q, q0)
% NDMATRIX(S, q, q0) returns the numerical derivative MAT of the given 
%   matrix S w.r.t. S at numerical value Q0
		
    r = size(q,1);
    [m, n] = size(S(q0));
    dS = repamat(zeros(m,n), r);
    
    
    for i = [1, r]
        % TODO: check what this is doing
        dS((i)) = substitute(S(SwitchAt(q0, q, i)), q(i), q0(i));
    end
end