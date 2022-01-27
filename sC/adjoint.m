function X = adjoint(Gst)
% ADJOINT(Gst) computes the adjoint transform X of matrix GST.
%   X is a 6x6 matrix that can be used congruent twists w.r.t. the
%   homogeneous transformation GST.

% TODO: add ADJOINT(rotm, disp)?
%   If 2 arguments are passed, they are considered as a matrix rotation
%   ROTM and a displacement vector DISP.
% nel caso di 2 input vengono assegnate la matrice di rotazione
% (rotm) e vett traslazione (disp) di Gst. (non serve mai)

R = Gst(1:3,1:3);
d = Gst(1:3,4);
X = [R,hat(d)*R;zeros(3,3),R];

% TODO: needed?
% for codegen with ADvar

% Rval = Gst.val(1:3,1:3);
% Rder = Gst.der(1:3,1:3);
% dval = Gst.val(1:3,4);
% dder = Gst.der(1:3,4);
% X = ADvar([Rval,hat(dval)*Rval;zeros(3,3),Rval], [Rder,hat(dder)*Rval+hat(dval)*Rder;zeros(3,3),Rder]);


end