function JacDot = bodyJacDerivative(jac, qdot)

%
% Js_dot = SpatialJac(Jac, qdot);
% Funzione che calcola la derivata dello Jacobiano body di un seriale

n = size(jac, 2);

JacDot = zeros(6, n);
adSum = zeros(6, 6); % adjoint cumulative sum initialization
for j = n:-1:1 % loop through the columns of the jacobian
    
    adSum = adSum - ad(jac(:, j)).*qdot(j); % adjoint of the previous columns
    JacDot(:, j) = adSum*jac(:, j);
    
end

end