function quat_v = QuatVectPart(mat)
    q = MatToQuat(mat);
    n = length(q);
    quat_v = q(2:n);
end