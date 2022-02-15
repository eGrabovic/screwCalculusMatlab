function JacDot = bodyJacDerivative(jac, qdot)
% BODYJACDERIVATIVE(Jac, qdot) computes the derivate JACDOT of the body 
%   Jacobian JAC of the robot w.r.t. the joint variables derivate QDOT.

    n = size(jac, 2);
    JacDot = zeros(6, n);

    % adjoint cumulative sum initialization
    adSum = zeros(6, 6); 
    
    % Loop through the columns of the jacobian
    for j = n:-1:1 
        % adjoint of the previous columns
        adSum = adSum - ad(jac(:, j)) .* qdot(j); 

        JacDot(:, j) = adSum * jac(:, j);   
    end

end