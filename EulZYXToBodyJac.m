function bjac = EulZYXToBodyJac(phi, theta, psi)
% EULZYXTOBODYJAC(alpha, beta, gamma) returns the Body Jacobian BJAC
%   corresponding to the ZYX Euler angles [PHI, THETA, PSI].
%   This is such that w_b = bjac * [phi, theta, psi], with w_b body fixed 
%   components. TODO: how to insert formula?

    bjac = [[           -sin(theta),         0,  1 ]         
		    [ cos(theta) * sin(psi),  cos(psi),  0 ]
		    [ cos(theta) * cos(psi), -sin(psi),  0 ]];
end