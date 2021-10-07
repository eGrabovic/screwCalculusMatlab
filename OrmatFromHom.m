function ormat = OrmatFromHom(mat)
%   Module[{nr, nc}]
    [nr, nc] = size(mat); 

    % Check to make sure that we were passed a square matrix *)
    if (nr ~= nc) || (nr < 4)
        print("Matrix has wrong dimensions");
        ormat = 0;
        return
    end

    % Extract the 3x3 upper left corner *)
    ormat = mat(1:3, 1:3);
end