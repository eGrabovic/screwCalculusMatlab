function [xi, theta] = LocalTranTwist(gStart, gEnd)
% Calculation of the error twist xi_err and error angle theta_err such that:
%     gStart.TwistExp(xi_err, theta_err) = gEnd

    assert(ismatrix(gstart), "First argument must be a matrix");
    assert(ismatrix(gEnd), "Second argument must be a matrix");
     
    gError = pinv(gStart) * gEnd;
    
    if all(gError == eye(4),'all')
        xi = zeros(1,6);
        theta = 0;
    else
        [xi,theta] = RigidTwist(gError);
    end
end