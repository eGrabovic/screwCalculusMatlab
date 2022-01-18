function Bmat = Inertia(DHtable, CGtable, Masslist, Tensortable)
    n = size(DHtable,1);
    
    Bmat = zeros(n);
    
    for k = [1, n]
    
        R0k  = RigidOrientation(DHFKine(DHtable, k));
        Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
        
        % attenzione: questo Jacobiano è relativo al baricentro
        
        Jvk =  Jk(1:3, :);
        Jok =  Jk(4:6, :);
        
        JvkT = Jvk';
        JokT = Jok';
        
        Bmat = Bmat + Masslist(k) * JvkT * Jvk + ...
                    JokT * R0k * Tensortable(k) * R0k' * Jok; 
    end
end