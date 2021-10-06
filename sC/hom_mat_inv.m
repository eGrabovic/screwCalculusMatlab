function GstInv = hom_mat_inv(Gst)
%
% Computes the inverse of a homogeneous transormation Gst into
% Gst^-1
%

Rst = Gst(1:3,1:3);
dst = Gst(1:3,4);
GstInv = [Rst.',-Rst.'*dst;0,0,0,1];

% for codegen with ADvar
% Rstval = Gst.val(1:3,1:3);
% Rstder = Gst.der(1:3,1:3);
% dstval = Gst.val(1:3,4);
% dstder = Gst.der(1:3,4);
% GstInv = ADvar([Rstval.',-Rstval.'*dstval;0,0,0,1], [Rstder.', -Rstder.'*dstval-Rstval.'*dstder; 0, 0, 0, 0]);
end