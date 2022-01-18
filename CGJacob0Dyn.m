function CGJacob = CGJacob0Dyn(DHtable, CGtable, k)
    R0k  = RigidOrientation(DHFKine(DHtable, k));
    pkck = CGtable(k);

    M = [[eye(3),        -Hat(R0k * pkck)];
        [zeros(3,3),    eye(3)]];

    DH_Jacob = DHJacob0Dyn(DHtable, k);
    CGJacob = M * DH_Jacob;
end