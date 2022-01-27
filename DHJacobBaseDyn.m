function jacobian = DHJacobBaseDyn(DHtable, Tb0, Tne, k)
% DHJACOBBASEDYN(DHtable, Tb0, Tne, k) computes the DH spatial JACOBIAN 
%   till the K-th joint with the needed empty columns stacked on the right,
%   considering also the initial offset TB0 and final offset TNE.

    n = size(DHtable, 2);
		     
    Jk  = [];    
    Hk  = DHFKine(DHtable, Tb0, Tne, k);
    
    for i = [1, k]
    
        type = DHtable(i, 5);
        zz = (RigidOrientation(DHFKine(DHtable, Tb0, Tne, i-1)));
        z = zz(:, 3);
        r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, Tb0, Tne, i-1)));
        
        if type == "P"        
            % Prismatic Joint
            jv = z;
            jo = zeros(1,3);
            j  = Join(jv, jo) ; 
        end
        
        if type == 'R'
            % Revolute Joint
            jv = Hat(z) * r;
            jo = z;
            j  = Join(jv, jo); 
        end

        Jk = [Jk, j];
    end
    
    jacobian = [Jk', zeros(6, n-k)];
end