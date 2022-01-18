function bjac = EulZYXToBodyJac(phi, theta, psi)
    bjac = [[           -sin(theta),         0,  1 ]         
		    [ cos(theta) * sin(psi),  cos(psi),  0 ]
		    [ cos(theta) * cos(psi), -sin(psi),  0 ]];
end