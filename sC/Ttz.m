function out = Ttz(z)
% TTX(x) computes the homogeneous matrix OUT describing a translation Z
%   along the z-axis.

    out = eye(4, class(z));

    out(3, 4) = z;

end