function ormat = OrmatFromHom(mat)
% ORMATFROMHOM(g) extracts rotation matrix ORMAT from the homogeneous
%   matrix MAT

% TODO: name changed for clarity; keep it?

    [nr, nc] = size(mat); 

    % Check to make sure that we were passed a square matrix
    if (nr ~= nc) || (nr < 4)
        print("Matrix has wrong dimensions");
        ormat = 0;
        return
    end

    % Extract the 3x3 upper left corner
    ormat = mat(1:3, 1:3);
end