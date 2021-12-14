function esp = TwistExpVector(xi, theta) 
    assert(isvector(xi), "This function only handles vectors");

    esp = TwistExp(HomogeneousToTwist(xi), theta); 
end