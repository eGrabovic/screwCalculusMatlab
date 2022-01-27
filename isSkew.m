function tf = isSkew(mat)
% ISSKEW(m) tests whether matrix M is a skew-symmetrix matrix and returns a
%   boolean accordingly

    if ~(ismatrix(mat)) || size(mat,1)~=size(mat,2)
        disp("Argument is not a square matrix.")
        tf = false;
    end

    tf = all(mat == mat');
end