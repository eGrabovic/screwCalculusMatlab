function jacobian = DHJacob0Dyn(DHtable, k)
% DHJACOB0DYN(DHtable, k) computes the DH spatial JACOBIAN till K-th joint 
%   with the needed empty columns stacked on the right, based on the
%   supplied Denavit-Hartemberg table DHTABLE.

    n = size(DHtable, 2);
    
    Jk  = [];     
    Hk = DHFKine(DHtable, k);
    
    for i = [1, k]
        
        type = DHtable(i, 5);
        zz = RigidOrientation(DHFKine(DHtable, i-1));
        z = zz(:,3);
        r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, i-1)));
        
        if type == "P" 
            % Prismatic Joint
            jv = z;
            jo = zeros(1,3);
            j  = Join(jv, jo) ;
        end
        
        if type == "R"
            % Revolute Joint
            jv = Hat(z) * r;
            jo = z;
            j  = Join(jv, jo); 
        end

        Jk = [Jk, j];
    end
    
    jacobian = [Jk', zeros(6, n-k)];
end