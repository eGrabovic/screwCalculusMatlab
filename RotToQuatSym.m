function quat = RotToQuatSym(mat)
% Warning: these are the reverse of conventional argument order:
% change order once all is consistent to (angle, nhat)
    
    tr = trace(mat);
    s = sqrt(tr + 1.0);
    q0 = s/2;
    s = 1/(2 * s);
    q1 = (mat(3,2) - mat(2,3)) * s;
    q2 = (mat(1,3) - mat(3,1)) * s;
    q3 = (mat(2,1) - mat(1,2)) * s;

    quat = [q0, q1, q2, q3];
end