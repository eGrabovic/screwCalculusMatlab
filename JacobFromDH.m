function jacobian = JacobFromDH(DHTable, upto)
% JACOBFROMDH builds the Jacobian of the robot related to the supplied
%     Denavit-Hartemberg DHTABLE. 
%
%     If UPTO contains a non-negative integer lower than the number of
%     joints inferrable from the table, the Jacobian built stops at the
%     UPTO-th joint.                                                    

    n = size(DHTable,1);
    J  = [];
    FK = DHFKine(DHTable);

    if exist('upto', 'var')
        if (upto > n) || (upto <= 0)
%             disp("Specified index out of buond. Using number of joints.");
            index = n;
        else
            index = upto;
        end
    else
        index = n;
    end
    
    % Build Jacobian joint-by-joint
    for i = [1:index]
        type = string(DHTable(i, 5));
        orient_mat = OrmatFromHom(DHFKine(DHTable, i-1));
        z =  orient_mat(:, 3);
        r = PosFromHom(FK) - (PosFromHom(DHFKine(DHTable, i-1)));
        
        if(type == 'P')
            % Prismatic Joint
            jac_v   = z; % Following DH convention, the joint axis is always 
                         % on the z-axis; that's the only case jac_v ~= 000
            jac_o   = [0; 0; 0];
            j       = [jac_v; jac_o];
        end
        
        if(type == 'R')
            % Revolute Joint
            jac_v   = Hat(z) * r;
            jac_o   = z;
            j       = [jac_v; jac_o]; 
        end
        
        J = [J, j];
    end
    
    jacobian = J;%[J, j];
end