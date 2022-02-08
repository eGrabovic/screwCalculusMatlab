function GstInv = rigidInverse(Gst)
% RIGIDINVERSE(Gst) computes the inverse GSTINV of a homogeneous matrix GST.
%   Allows ADvar class Gst argument.

    % Manage ADvar argument
    if isa(Gst, 'ADvar')
        Der = rigidInverse(Gst.der);
        GstInv = ADvar(rigidInverse(Gst.val), Der);
        return 
    end


    assert(all(size(mat) == [4 4]), ...
        "Provided matrix must be an ADvar instance or a 4x4 homogeneous matrix");

    Rst = Gst(1:3, 1:3);
    dst = Gst(1:3, 4);
    GstInv = [Rst.', -Rst.'*dst; 0, 0, 0, Gst(4, 4)];
end