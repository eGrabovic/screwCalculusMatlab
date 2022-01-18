function Gvec = Gravitational(DHtable, CGtable, Masslist, g0)
    n = size(DHtable,1);
	
	Gmat = zeros(1, n);
    Gvec = Gmat(1);
	
    for k = [1, n]
		
		%R0k  = RigidOrientation(DHFKine(DHtable, k));
	    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
	    
	    % attenzione: questo Jacobiano è relativo al baricentro
	    
	    Jvk =  Jk(1:3, :);
	    % Jok =  Jk(4:6, :);
	    
	    JvkT = Jvk';
	    % JokT = Jok';
	    
	    Gvec = Gvec - Masslist(k) * JvkT * g0; 
    end
end