function JacDot = spatialJacDerivative(jac, qdot)
% SPATIALJACDERIVATIVE(jac, qdot) computes the derivate JACDOT of the 
%   spatial Jacobian JAC of a robot w.r.t. the joint variables derivate QDOT.

    n = size(jac, 2);
    JacDot = zeros(6, n);

    % adjoint cumulative sum initialization
    adSum = zeros(6, 6); 

    % Loop through the columns of the jacobian
    for j = 2:n 
        % adjoint of the previous columns
        adSum = adSum + ad(jac(:, j-1)).*qdot(j-1); 

        JacDot(:, j) = adSum*jac(:, j);
    end

end