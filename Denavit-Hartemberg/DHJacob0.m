function jacobian = DHJacob0(DHTable, upto)
% DHJACOB0(DHTable, upto) builds the JACOBIAN of the robot related to the 
%   supplied Denavit-Hartemberg DHTABLE. 
%
%   If UPTO contains a non-negative integer lower than the number of
%   joints inferrable from the table, the Jacobian built stops at the
%   UPTO-th joint.                                                    
	
    n = size(DHTable,1);
    J  = [];
    FK = DHFKine(DHTable);

    if exist('upto', 'var')
        assert((upto < n) && (upto >= 0), ...
            "Specified index must be non-negative and lower than " + ...
            "the number of joints")

        if upto == 0 % Remaining edge case
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
    
    jacobian = J;
    
    % If the Jacobian is referred to a non-endeffector joint, pad the
    % matrix to match the dimensions of the body
    if index ~= n
        jacobian = [jacobian, repmat(0, size(jacobian, 1), n-index)];
    end
end