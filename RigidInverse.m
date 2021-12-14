function rig_inv = RigidInverse(mat)
% Calculate the inverse rigid transformation

    assert(ismatrix(mat), "Provided argument must be a matrix");
    assert(all(size(mat) == [4 4]), ...
        "Provided matrix must be a 4x4 homogeneous matrix");

    Rot = RigidOrientation(mat); 
    pos = RigidPosition(mat);

    rig_inv = RPToHomogeneous(Rot', -Rot' * pos);
end