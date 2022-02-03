function xi_dim = xidim(xi)
% XIDIM returns the space dimension XI_DIM of a twist vector XI.
%   Input
%       xi:     twist vector
%   Output
%       xi_dim: space dimension of xi
%
%   Example
%       xidim(zeros(6,1)) --> 3      (meaning xi is in R3)
%       xidim(zeros(3,1)) --> 2

    assert(isvector(xi), "xi is not a vector");

    l = max(size(xi));

    % Check the dimensions of the vector to make sure everything is okay
    n = (sqrt(1 + 8*l) - 1)/2;

    assert(~isinf(n) && (floor(n) == n), "The input xi has wrong dimensions");

    xi_dim = n;
end