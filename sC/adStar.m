function out = adStar(twist)
% ADSTAR(twist) computes the transpose X of the adjoint trasformation of
%   the matrix corresponding to the twist array TWIST.

    w = twist(4:6);
    v = twist(1:3);
    what = hat(w);

    out = [[what,   zeros(3, 3)];
           [hat(v), what]];

end