function Cmat = InertiaToCoriolisNotChristoffel(M, q, qp)
% INERTIATOCORIOLISCHRISTOFFERL(M, q, qp) computes the Coriolis matrix CMAT
%   given:
%   - the inertia matrix B
%   - a list of the joint variables Q
%   - the list of joint variable derivatives QP

    n = size(M,1);
    Cmat = zeros(n,n);
    
    
    for i = [1, n]
        for j = [1, n]
            for k = [1, n]
            Cmat(i,j) = Cmat(i,j) + ...
                (diff(M(i, j), q(k)) - (1/2) * diff(M(j, k), q(i)))* qp(k);
            end
        end
    end
end