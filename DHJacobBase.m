function jacobian = DHJacobBase(DHtable, Tb0)
    R    = RigidOrientation(Tb0);
    g    = RPToHomogeneous(R, zeros(1,3));
    Adg  = RigidAdjoint(g);
    J    = DHJacob0(DHtable);
    jacobian   = Adg * J;
end