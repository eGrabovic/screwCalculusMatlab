function ang_twist = xitow(xi)
% XITOW(xi) gives the angular part ANG_TWIST of twist XI.

% Extract the angular portion of a twist (private)

    assert(isvector(xi), "Input xi must be a vector");

    n = xidim(xi);

    % Make sure that the vector had a reasonable length
    assert(n > 0, "Input array cannot be empty");

    % Extract the angular portion of the twist
    ang_twist = xi(end - (n * (n-1) / 2 - 1):end);
end