function pos = RigidPosition(mat)
% RIGIDPOSITION(mat) extracts the position vector POS from matrix MAT
  
    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2) || size(mat,1) < 4
        disp("Argument is not a 4x4 homogeneous matrix")
        return
    end

    pos = mat(1:3, 4);
end