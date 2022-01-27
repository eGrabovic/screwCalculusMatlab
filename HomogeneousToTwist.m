function twist = HomogeneousToTwist(A)
% HOMOGENEOUSTOTWIST(A) converts A from a 4x4 matrix to a 6 element vector
%   TWIST.
%! This only works in dimensions 2 and 3 for now!

    % Check to make sure that our input makes sense
    assert(ismatrix(A), "A is not a matrix");
    [nr, nc] = size(A); 
    assert(nr == nc, "A is not square");

    % Make sure that we have a twist and not a rigid motion
    assert(A(nr,nc) == 0, ...
        "The matrix represents a rigid motion, not a twist");

    % Extract the skew part and the vector part and make a vector
    twist = [A(1:nr, nr-1:1), SkewToAxis(A(1,nr-1))];
end