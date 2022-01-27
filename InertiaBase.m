function Bmat = InertiaBase(DHtable, CGtable, Masslist, Tensortable, Tb0, Tne)
% INERTIABASE(DHtable, CGtable, Masslist, Tensortable, Tb0, Tne) returns 
%   matrix BMAT [B(q)] for the supplied: 
%   - DHTABLE (Denavit-Hartemberg table)
%   - GCTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - TENSORTABLE(table with link inertia tensors)
%   - [TB0, TNE] (initial [B to S0] and final [Sn to E] offset transformations)

    n = size(DHtable, 1);
    
    Bmat = zeros(n);
    
    for k = [1, n]
    
        R0k  = RigidOrientation(DHFKine(DHtable, Tb0, Tne, k));
        Jk  =  CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k);
        
        % WARNING: this Jacobian refers to the barycenter
        
        Jvk =  Jk(1:3, :);
        Jok =  Jk(4:6, :);
        
        JvkT = Jvk';
        JokT = Jok';
        
        Bmat = Bmat + Masslist(k) * JvkT * Jvk + ...
                    JokT * R0k * Tensortable(k) * R0k' * Jok; 
    end
end