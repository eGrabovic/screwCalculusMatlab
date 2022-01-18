function Gvec = GravitationalBase(DHtable, CGtable, Masslist, g0, Tb0, Tne)
	n = size(DHtable,1);
		
    Gmat = zeros(1, n);
    Gvec = Gmat(1);
		
    for k = [1, n]
			
		% R0k  = RigidOrientation( DHFKine(DHtable, Tb0, Tne, k));
	    Jk  =  CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k);
	    
	    % attenzione: questo Jacobiano è relativo al baricentro
	    
	    Jvk =  Jk(1:3, :);
	    % Jok =  Jk(4:6, :);
	    
	    JvkT = Jvk';
	    % JokT = Jok';
	    
	    Gvec = Gvec - Masslist(k) * JvkT * g0; 
    end
end