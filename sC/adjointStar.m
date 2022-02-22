function X = adjointStar(Gst)
% ADJOINTSTAR(Gst) computes the transpose X of the adjoint trasformation of
%   the matrix Gst.

    R = Gst(1:3, 1:3);
    d = Gst(1:3, 4);

    X = [[R,            zeros(3, 3)];
         [hat(d) * R,   R]];

end