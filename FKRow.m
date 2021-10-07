function hom_FK = FKRow(table_row) %DH
    a = table_row(1);
    alpha = table_row(2);
    d = table_row(3);
    theta = table_row(4);
%     type = table_row(4);
    
    % TODO: fix homogeneous; all to hom
    Tz = HomZ(0, [0, 0, d]);
    Rz = HomZ(theta, [0, 0, 0]); %Rotz(theta);
    Tx = HomX(0, [a, 0, 0]);
    Rx = HomX(alpha, [0, 0, 0]); %RotX(alpha);
    
    hom_FK = Tz*Rz*Tx*Rx;
end