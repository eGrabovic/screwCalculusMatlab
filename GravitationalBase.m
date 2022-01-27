function Gvec = GravitationalBase(DHtable, CGtable, Masslist, g0, Tb0, Tne)
% GRAVITATIONALBASE(DHtable, CGtable, Masslist, gb, Tb0, Tne) computes the 
%   gravitational vector GVEC for the supplied: 
%   - DHTABLE (Denavit-Hartemberg table)
%   - CGTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - GB (componenents of gravity in {B})
%   - [Tb0, Tne] initial [B to S0] and final [Sn to E] offset transformations

	n = size(DHtable,1);
		
    Gmat = zeros(1, n);
    Gvec = Gmat(1);
		
    for k = [1, n]
			
		% R0k  = RigidOrientation( DHFKine(DHtable, Tb0, Tne, k));
	    Jk  =  CGJacobBaseDyn(DHtable, CGtable, Tb0, Tne, k);
	    
        % WARNING: thia Jacobian is w.r.t. the barycenter
	    Jvk =  Jk(1:3, :);
	    % Jok =  Jk(4:6, :);
	    
	    JvkT = Jvk';
	    % JokT = Jok';
	    
	    Gvec = Gvec - Masslist(k) * JvkT * g0; 
    end
end