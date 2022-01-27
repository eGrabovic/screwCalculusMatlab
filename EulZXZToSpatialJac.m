function sjac = EulZXZToSpatialJac(phi, theta, psi)
% EULZXZTOSPATIALJAC(alpha, beta, gamma) returns the Spatial Jacobian SJAC
%   corresponding to the ZXZ Euler angles [PHI, THETA, PSI].
%   This is such that w_s = sjac * [phi, theta, psi], with w_s spatial 
%   components. TODO: how to insert formula?

    sjac = [[0, cos(phi),  sin(phi) * sin(theta)]
           [0, sin(phi), -cos(phi) * sin(theta)]
           [1,        0,            cos(theta)]];
end