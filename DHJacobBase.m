function jacobian = DHJacobBase(DHtable, Tb0)
% DHJACOBBASE(DHtable, Tb0) computes the DH spatial JACOBIAN for the 
%   supplied DH table DHTABLE w.r.t. a base frame Sb (both origin and 
%   orientation are considered), where the displacement from Sb to S0 is 
%   expressed by the homogeneous matrix TB0.

    R    = RigidOrientation(Tb0);
    g    = RPToHomogeneous(R, zeros(1,3));
    Adg  = RigidAdjoint(g);
    J    = DHJacob0(DHtable);
    jacobian   = Adg * J;
end