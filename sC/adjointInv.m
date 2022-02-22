function X = adjointInv(Gst)
% ADJOINTINV(Gst) computes the inverse X of the adjoint transformation of 
%   the matrix GST.

    R = Gst(1:3, 1:3).';
    d = Gst(1:3, 4);

    X = [[R,            -R * hat(d)];
         [zeros(3, 3),  R]];

end