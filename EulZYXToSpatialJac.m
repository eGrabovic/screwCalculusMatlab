function sjac = EulZYXToSpatialJac(phi, theta, psi)
% EULZYXTOSPATIALJAC(alpha, beta, gamma) returns the Spatial Jacobian SJAC
%   corresponding to the ZYX Euler angles [PHI, THETA, PSI].
%   This is such that w_s = sjac * [phi, theta, psi], with w_s spatial 
%   components. TODO: how to insert formula?

    sjac = [[ 0, -sin(phi), cos(theta) * cos(phi)]
   		    [ 0,  cos(phi), cos(theta) * sin(phi)]
		    [ 1,        0,            -sin(theta)]];
end