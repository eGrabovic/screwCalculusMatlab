function sjac = EulZYXToSpatialJac(phi, theta, psi)
    sjac = [[ 0, -sin(phi), cos(theta) * cos(phi)]
   		    [ 0,  cos(phi), cos(theta) * sin(phi)]
		    [ 1,        0,            -sin(theta)]];
end