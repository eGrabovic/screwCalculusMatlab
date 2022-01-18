function CGJacob = CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k)
    R0k  = RigidOrientation(DHFKine(DHtable, Tb0, Tne, k));
    pkck = CGtable(k);

    M = [[eye(3),     -Hat(R0k * pkck)];
        [zeros(3,3),  eye(3)]];
    
    DH_Jacob = DHJacobBaseDyn(DHtable, Tb0, Tne, k);
    CGJacob = M * DH_Jacob;
end