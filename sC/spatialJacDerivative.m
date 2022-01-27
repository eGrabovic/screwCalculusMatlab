function JacDot = spatialJacDerivative(jac, qdot)
% SPATIALJACDERIVATIVE(jac, qdot) computes the derivate JACDOT of the 
%   spatial Jacobian of a robot.

n = size(jac, 2);

JacDot = zeros(6, n);
adSum = zeros(6, 6); % adjoint cumulative sum initialization
for j = 2:n % loop through the columns of the jacobian
    
    adSum = adSum + ad(jac(:, j-1)).*qdot(j-1); % adjoint of the previous columns
    JacDot(:, j) = adSum*jac(:, j);
    
end

end