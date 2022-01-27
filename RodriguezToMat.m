function mat = RodriguezToMat(gamma1, gamma2, gamma3)
% RODRIGUEZTOMAT(gamma1, gamma2, gamma3) returns the rotation matrix MAT
%   corresponding to the supplied Rodriguez parameters GAMMA1, GAMMA2,
%   GAMMA3.
%   Recall that 
%   {g1, g2, g3} = r tan(theta/2)
%   TODO: how to add formulas?

% Rodriguez parameters gamma = r * tan(theta/2)
    gamma = [gamma1, gamma2, gamma3];
    hatgamma = Hat(gamma);
    modulusgammasquared = gamma * gamma;
    
    mat = eye(3) + 2 / (1 + modulusgammasquared) * ...
                       (hatgamma + hatgamma * hatgamma);
end