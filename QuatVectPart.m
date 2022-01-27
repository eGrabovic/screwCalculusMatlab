function quat_v = QuatVectPart(mat)
% QUATVECTPART(MAT) returns the 3-vector part QUAT_V of quaternion (q0, q) 
%   associated with the rotation matrix MAT

    q = MatToQuat(mat);
    n = length(q);
    quat_v = q(2:n);
end