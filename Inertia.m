function Bmat = Inertia(DHtable, CGtable, Masslist, Tensortable)
% INERTIA(DHtable, CGtable, Masslist, Tensortable) returns matrix BMAT
%   [B(q)] for the supplied: 
%   - DHTABLE (Denavit-Hartemberg table)
%   - GCTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - TENSORTABLE (table with link inertia tensors)

    n = size(DHtable,1);
    
    Bmat = zeros(n);
    
    for k = [1, n]
    
        R0k  = RigidOrientation(DHFKine(DHtable, k));
        Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
        
        % WARNING: this Jacobian refers to the barycenter
        
        Jvk =  Jk(1:3, :);
        Jok =  Jk(4:6, :);
        
        JvkT = Jvk';
        JokT = Jok';
        
        Bmat = Bmat + Masslist(k) * JvkT * Jvk + ...
                    JokT * R0k * Tensortable(k) * R0k' * Jok; 
    end
end