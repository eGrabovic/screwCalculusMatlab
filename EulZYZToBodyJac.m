function bjac = EulZYZToBodyJac(phi, theta, psi)
% EULZYZTOBODYJAC(alpha, beta, gamma) returns the Body Jacobian BJAC 
%   corresponding to the ZYZ Euler angles [PHI, THETA, PSI].
%   This is such that w_b = bjac * [phi, theta, psi], with w_b body fixed
%   components. TODO: how to insert formula?

		bjac = [[ -cos(psi) * sin(theta), sin(psi), 0]
		        [  sin(psi) * sin(theta), cos(psi), 0]
		        [             cos(theta),        0, 1]];
end