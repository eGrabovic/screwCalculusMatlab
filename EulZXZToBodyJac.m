function bjac = EulZXZToBodyJac(phi, theta, psi)
% EULTOZXZBODYJAC(alpha, beta, gamma) returns the Body Jacobian BJAC
%   corresponding to the ZXZ Euler angles [PHI, THETA, PSI].
%   This is such that w_b = bjac * [phi, theta, psi], with w_b body fixed
%   components. TODO: how to insert formula?

    bjac = [[ sin(theta) * sin(psi),  cos(psi), 0 ]
		    [ sin(theta) * cos(psi), -sin(psi), 0 ]
		    [            cos(theta),         0, 1 ]];
end