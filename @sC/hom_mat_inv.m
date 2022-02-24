function GstInv = hom_mat_inv(Gst)
% HOM_MAT_INV(Gst) computes the inverse GSTINV of a homogeneous matrix GST.

    if all(size(Gst) ~= [4,4])
        error('Homogenous transformation matrix must have 4x4 dimension')
    end

    Rst = Gst(1:3, 1:3);
    dst = Gst(1:3, 4);

    GstInv = [[Rst',        -Rst' * dst];
              [0, 0, 0,     1]];
end