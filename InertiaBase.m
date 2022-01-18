function Bmat = InertiaBase(DHtable, CGtable, Masslist, Tensortable, Tb0, Tne)
    n = size(DHtable, 1);
    
    Bmat = zeros(n);
    
    for k = [1, n]
    
        R0k  = RigidOrientation(DHFKine(DHtable, Tb0, Tne, k));
        Jk  =  CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k);
        
        % attenzione: questo Jacobiano è relativo al baricentro
        
        Jvk =  Jk(1:3, :);
        Jok =  Jk(4:6, :);
        
        JvkT = Jvk';
        JokT = Jok';
        
        Bmat = Bmat + Masslist(k) * JvkT * Jvk + ...
                    JokT * R0k * Tensortable(k) * R0k' * Jok; 
    end
end