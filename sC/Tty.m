function out = Tty(y)
% TTX(x) computes the homogeneous matrix OUT describing a translation Y
%   along the y-axis.

    out = eye(4, class(y));

    out(2, 4) = y;

end