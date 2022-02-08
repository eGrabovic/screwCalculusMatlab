function GstInv = hom_mat_inv(Gst)
% HOM_MAT_INV(Gst) computes the inverse GSTINV of a homogeneous matrix GST.

    assert(all(size(Gst) == [4 4]), ...
        "Provided matrix must be a 4x4 homogeneous matrix");

    Rst = Gst(1:3, 1:3);
    dst = Gst(1:3, 4);

    GstInv = [Rst.', -Rst.' * dst;
              0, 0, 0, 1];
end