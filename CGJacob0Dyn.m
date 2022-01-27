function CGJacob = CGJacob0Dyn(DHtable, CGtable, k)
% CGJACOB0DYN(DHtable, CGtable, k) computes the CG (linear velocity of the 
%   CG) spatial Jacobian till the K-th joint with the needed empty columns 
%   stacked on the right
%
%   Besides the DH table DHTABLE, a table with the pici (CGTABLE) must be 
%   supplied.

    R0k  = RigidOrientation(DHFKine(DHtable, k));
    pkck = CGtable(k);

    M = [[eye(3),        -Hat(R0k * pkck)];
        [zeros(3,3),    eye(3)]];

    DH_Jacob = DHJacob0Dyn(DHtable, k);
    CGJacob = M * DH_Jacob;
end