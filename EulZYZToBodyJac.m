function bjac = EulZYZToBodyJac(phi, theta, psi)
		bjac = [[ -cos(psi) * sin(theta), sin(psi), 0]
		        [  sin(psi) * sin(theta), cos(psi), 0]
		        [             cos(theta),        0, 1]];
end