function [xi, theta] = LocalTranTwist(gStart, gEnd)
% LOCALTRANTWIST(gStart, gEnd) returns the twist XI and angle THETA such 
%   that: 
%   gStart * TwistExp(xi,theta) = gEnd
%   The xi coordinates are expressed in the local (gStart) reference frame.
% TODO: how to add formulas?

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