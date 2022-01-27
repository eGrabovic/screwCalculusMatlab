function rot = RigidOrientation(mat)
% RIGIORIENTATION(mat) extracts rotation matrix ROT from matrix MAT

    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2) || size(mat,1) < 4
        disp("Argument is not a 4x4 homogeneous matrix")
        return
    end

    rot = mat(1:3, 1:3);
end