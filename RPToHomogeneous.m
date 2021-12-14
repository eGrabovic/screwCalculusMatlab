function hom = RPToHomogeneous(Rot, pos)
% Convert a rotation + a translation to a homogeneous matrix *)

    % Check to make sure the dimensions of the args make sense *)
    assert(isvector(pos), "pos is not a vector");
    assert(ismatrix(Rot), "Rot is not a matrix");

    n = max(size(pos));

    assert(n == 3, ...
        "Array size is not acceptable; it must consist of 3 elements");
    assert(all(sizes(Rot) == [n, n]), ...
        "Matrix and vector sizes are not consistent");
	
    % Now put everything together into a homogeneous transformation *)
    hom = [[Rot,       pos];
           [zeros(1,3),  1]];
end