function pos = PosFromHom(mat)
% RIGIDPOSITION(g) extracts position vector POS from the matrix MAT
%   TODO: changed name for clarity; keep it?
    [nr, nc] = size(mat);
    
    % Check to make sure that we were passed a square matrix
    if (nr ~= nc) || (nr < 4)
        print("Matrix has wrong dimensions");
        pos = 0;
        return
    end

    % Extract the upper left column
    pos = mat(1:3, nc);
end