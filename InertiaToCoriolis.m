function Cmat = InertiaToCoriolis(M, q, qp)
% INERTIATOCORIOLIS(M, q, qp) computes the Coriolis matrix CMAT given:
%   - the inertia matrix B
%   - the joint variables Q
%   - the joint variable derivatives QP

    n = size(M,1);

    % Brute force calculation
    Cmat = zeros(n);

    for i = [1, n]
      for j = [1, n]
        for k = [1, n]
          Cmat(i,j) = Cmat(i,j) + (1/2) * qp(k) * diff(M(i,j), q(k)) + ...
                        diff(M(i,k), q(j)) - diff(M(j,k), q(i));
        end
      end
    end
end