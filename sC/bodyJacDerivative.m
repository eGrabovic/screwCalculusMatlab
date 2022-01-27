function JacDot = bodyJacDerivative(jac, qdot)
% TODO: why does this refer to spatialjac?
% SPATIALJAC(Jac, qdot) computes the derivate of the body Jacobian JAC of
%   the robot w.r.t. the joint variables derivate QDOT.

n = size(jac, 2);

JacDot = zeros(6, n);
adSum = zeros(6, 6); % adjoint cumulative sum initialization
for j = n:-1:1 % loop through the columns of the jacobian
    
    adSum = adSum - ad(jac(:, j)).*qdot(j); % adjoint of the previous columns
    JacDot(:, j) = adSum*jac(:, j);
    
end

end