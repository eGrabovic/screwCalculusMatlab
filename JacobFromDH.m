function jacobian = JacobFromDH(DHTable)
    n = size(DHTable,1);
    J  = [];
    FK = ForwardKinematics(DHTable);
    
    for i = [1:n]
        type = string(DHTable(i, 5)); % DHTable(i, 5);
        orient_mat = OrmatFromHom(ForwardKinematics(DHTable, i-1));
        z =  orient_mat(:, 3);
        r = PosFromHom(FK) - (PosFromHom(ForwardKinematics(DHTable, i-1)));
        
        if(type == 'P')
            % Prismatic Joint
            jac_v = z; % seguendo DH asse sempre su z: unico caso non 000
            jac_o = [0; 0; 0];
            j  = [jac_v; jac_o];
        end
        
        if(type == 'R')
            % Revolute Joint
            jac_v = Hat(z) * r;
            jac_o = z;
            j  = [jac_v; jac_o]; 
        end
        
        J = [J, j];
    end
    
    jacobian = J;%[J, j];
end