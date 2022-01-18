function mat = RodriguezToMat(gamma1, gamma2, gamma3)
% Rodriguez parameters gamma = r * tan(theta/2)
    gamma = [gamma1, gamma2, gamma3];
    hatgamma = Hat(gamma);
    modulusgammasquared = gamma * gamma;
    
    mat = eye(3) + 2 / (1 + modulusgammasquared) * ...
                       (hatgamma + hatgamma * hatgamma);
end