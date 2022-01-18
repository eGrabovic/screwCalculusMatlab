function sjac = EulZXZToSpatialJac(phi, theta, psi)
    sjac = [[0, cos(phi),  sin(phi) * sin(theta)]
           [0, sin(phi), -cos(phi) * sin(theta)]
           [1,        0,            cos(theta)]];
end