function lin_twist = xitov(xi)
% Extract the linear portion of a twist (private)

    assert(isvector(xi), "Provided xi must be a vector");
    n = xidim(xi);

    % Make sure that the vector had a reasonable length
    assert(n > 0, "Input has wrong dimensions");

    % Extract the linear portion of the twist
    lin_twist = xi(1:n);
end