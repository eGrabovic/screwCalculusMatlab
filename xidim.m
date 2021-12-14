function xi_dim = xidim(xi)
% Figure out the space dimension of a twist (private)

    assert(isvector(xi), "xi is not a vector");

    l = max(size(xi));

    % Check the dimensions of the vector to make sure everything is okay
    n = (sqrt(1 + 8*l) - 1)/2;

    assert(~isinf(n) && (floor(n) == n), "The input xi has wrong dimensions");

    xi_dim = n;
end