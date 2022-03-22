function GstInv = hom_mat_inv2D(screwObj, Gst)
% HOM_MAT_INV2D(Gst) computes the inverse GSTINV of a homogeneous matrix 
%   GST in R2.

    Rst = Gst(1:2, 1:2);
    dst = Gst(1:2, 3);

    GstInv = [[Rst.',   -Rst.' * dst];
              [0, 0,     1]];
    
end