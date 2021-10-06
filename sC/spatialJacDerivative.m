function JacDot = spatialJacDerivative(jac, qdot)

%
% Js_dot = SpatialJac(Jac, qdot);
% Funzione che calcola la derivata dello Jacobiano spatial di un seriale

n = size(jac, 2);

JacDot = zeros(6, n);
adSum = zeros(6, 6); % adjoint cumulative sum initialization
for j = 2:n % loop through the columns of the jacobian
    
    adSum = adSum + ad(jac(:, j-1)).*qdot(j-1); % adjoint of the previous columns
    JacDot(:, j) = adSum*jac(:, j);
    
end

end