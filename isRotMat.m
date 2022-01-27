function tf = isRotMat(mat)
% ISROTMAT(m) tests whether matrix M is a rotation matrix and returns a
%   boolean accordingly
% TODO: improve

    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2)
        disp("Argument is not a square matrix")
        tf = false;
    end

    nr = size(mat,1);

    tf = all(mat' * mat == eye(nr));
end