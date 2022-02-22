function X = adjoint(Gst)
% ADJOINT(Gst) computes the adjoint transformation X of matrix GST.
%   X is a 6x6 matrix that can be used as a congruent twists w.r.t. the
%   homogeneous transformation GST.

    R = Gst(1:3, 1:3);
    d = Gst(1:3, 4);

    X = [[R,            hat(d) * R];
        [zeros(3, 3),   R]];

end