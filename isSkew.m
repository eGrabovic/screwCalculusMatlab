function tf = isSkew(mat)
    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2)
        disp("Argument is not a square matrix.")
        tf = false;
    end

    tf = all(mat == mat');
end