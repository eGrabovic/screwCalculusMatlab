function jac = RodriguezToSpatialJac(gamma1, gamma2, gamma3)
% RODRIGUEZTOSPATIALJAC(gamma1, gamma2, gamma3) returns the Spatial 
%   Jacobian JAC corresponding to the supplied Rodriguez parameters GAMMA1,
%   GAMMA2, GAMMA3.
%   This is such that 
%   w_s = jac * dot{g1, g2, g3}
%   with w_s spatial components
%   TODO: how to add formulas?

    gamma = [gamma1, gamma2, gamma3];
    modulusgammasquared = gamma * gamma;

    jac = 2 / (1 + modulusgammasquared) *...
              [[1,          -gamma3,    gamma2];
			   [gamma3,     1,          -gamma1];
			   [-gamma2,    gamma1,     1]];
end