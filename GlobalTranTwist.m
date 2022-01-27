function [xi, theta] = GlobalTranTwist(gStart, gEnd)
% GLOBALTRANTWIST(gStart, gEnd) returns the twist error XI and angle error 
%   THETA such that: TwistExp(xi,theta) * gStart = gEnd.
%   TODO: how to insert formula?
%   The xi coordinate are expressed in the global (identity) reference frame.

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
