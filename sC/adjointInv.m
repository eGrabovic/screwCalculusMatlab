function X = adjointInv(Gst)
%
% adjoint(Gst)
% trasformata aggiunta di Gst. Matrice 6x6 che può trasformare
% twist congruenti in base alla trasformazione omogenea di Gst.
% adjoint(rotm,disp)
% nel caso di 2 input vengono assegnate la matrice di rotazione
% (rotm) e vett traslazione (disp) di Gst. (non serve mai)
%

R = Gst(1:3,1:3).';
d = Gst(1:3,4);
X = [R,-R*hat(d);zeros(3,3),R];

% for codegen with ADvar

% Rval = Gst.val(1:3,1:3);
% Rder = Gst.der(1:3,1:3);
% dval = Gst.val(1:3,4);
% dder = Gst.der(1:3,4);
% X = ADvar([Rval,hat(dval)*Rval;zeros(3,3),Rval], [Rder,hat(dder)*Rval+hat(dval)*Rder;zeros(3,3),Rder]);


end