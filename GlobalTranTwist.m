function [xi, theta] = GlobalTranTwist(gStart, gEnd)
% Calculation of the error twist xi_err and error angle theta_err such that:
% TwistExp(xi_err, theta_err).gStart = gEnd

    assert(ismatrix(gStart), "First argument must be a matrix");
    assert(ismatrix(gEnd), "First argument must be a matrix");
     
     gError = pinv(gStart) * gEnd;
     
     if all(gError == eye(4), 'all')
       xi = zeros(1,6);
       theta = 0;
     else
         [xi, theta] = RigidTwist(gError);
         
         Adg = RigidAdjoint(gStart);
         xi = Adg * xi;
     end
end
