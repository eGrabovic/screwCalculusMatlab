function bjac = EulZXZToBodyJac(phi, theta, psi)
    bjac = [[ sin(theta) * sin(psi),  cos(psi), 0 ]
		    [ sin(theta) * cos(psi), -sin(psi), 0 ]
		    [            cos(theta),         0, 1 ]];
end