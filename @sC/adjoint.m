function X = adjoint(rotm,disp)
    %
    % adjoint(Gst)
    % trasformata aggiunta di Gst. Matrice 6x6 che può trasformare
    % twist congruenti in base alla trasformazione omogenea di Gst.
    % adjoint(rotm,disp)
    % nel caso di 2 input vengono assegnate la matrice di rotazione
    % (rotm) e vett traslazione (disp) di Gst.
    %
    if nargin == 1
        if all(size(rotm) ~= [4,4])
            error('size of Gst homogeneous matrix requires to be 4x4')
        end
        R = rotm(1:3,1:3);
        d = rotm(1:3,4);
        X = [R,sC.hat(d)*R;zeros(3,3),R];
    else
        if all(size(rotm) ~= [3,3]) || all(size(disp) ~= [3,1])
            error('rotation matrix requires to be a 3x3 matrix and the displacement vector requires to be a 3x1 matrix')
        end
        R = rotm;
        d = disp;
        X = [R sC.hat(d)*R;zeros(3,3) R];
    end
end