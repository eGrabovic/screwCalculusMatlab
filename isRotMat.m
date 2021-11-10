function tf = isRotMat(mat)
    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2)
        disp("Argument is not a square matrix")
        tf = false;
    end

    nr = size(mat,1);

    tf = all(mat' * mat == eye(nr));
end