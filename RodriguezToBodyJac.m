function jac = RodriguezToBodyJac(gamma1, gamma2, gamma3)
		gamma = [gamma1, gamma2, gamma3];
		modulusgammasquared = gamma * gamma;

		jac = 2 / (1 + modulusgammasquared) *...
                  [[1,          gamma3,    -gamma2];
				   [-gamma3,    1,          gamma1];
				   [gamma2,    -gamma1,     1]];
end