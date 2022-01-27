function CGJacob = CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k)
% CGJACOBBASEDYN(DHtable, CGtable, Tb0, Tne, k) computes the CG (linear 
%   velocity of the CG) spatial Jacobian till the K-th joint with the 
%   needed empty columns stacked on the right.
%
%   Besides the DH table DHTABLE, a table with the pici (CGTABLE) and 
%   the homogeneous matrices TB0, TNE (initial [B to S0] and final [Sn to 
%   E] offset transformations) must be supplied.

    R0k  = RigidOrientation(DHFKine(DHtable, Tb0, Tne, k));
    pkck = CGtable(k);

    M = [[eye(3),     -Hat(R0k * pkck)];
        [zeros(3,3),  eye(3)]];
    
    DH_Jacob = DHJacobBaseDyn(DHtable, Tb0, Tne, k);
    CGJacob = M * DH_Jacob;
end