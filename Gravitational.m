function Gvec = Gravitational(DHtable, CGtable, Masslist, g0)
% GRAVITATIONAL(DHtable, CGtable, Masslist, g0) computes the gravitational 
%   vector GVEC for the supplied: 
%   - DHTABLE (Denavit-Hartemberg table)
%   - CGTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - G0 (componenents of gravity in {0})

    n = size(DHtable,1);
	
	Gmat = zeros(1, n);
    Gvec = Gmat(1);
	
    for k = [1, n]
		
		%R0k  = RigidOrientation(DHFKine(DHtable, k));
	    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
	    
	    % WARNING: this Jacobian is w.r.t the barycenter
	    Jvk =  Jk(1:3, :);
	    % Jok =  Jk(4:6, :);
	    
	    JvkT = Jvk';
	    % JokT = Jok';
	    
	    Gvec = Gvec - Masslist(k) * JvkT * g0; 
    end
end