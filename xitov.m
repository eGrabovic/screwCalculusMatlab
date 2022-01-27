function lin_twist = xitov(xi)
% XITOV(xi) gives the translational part LIN_TWIST of twist XI.

% Extract the linear portion of a twist (private)

    assert(isvector(xi), "Provided xi must be a vector");
    n = xidim(xi);

    % Make sure that the vector had a reasonable length
    assert(n > 0, "Input has wrong dimensions");

    % Extract the linear portion of the twist
    lin_twist = xi(1:n);
end