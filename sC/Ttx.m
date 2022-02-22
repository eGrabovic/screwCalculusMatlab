function out = Ttx(x)
% TTX(x) computes the homogeneous matrix OUT describing a translation X
%   along the x-axis.

    out = eye(4, class(x));

    out(1, 4) = x;

end