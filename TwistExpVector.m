function esp = TwistExpVector(xi, theta) 
% TWISTEXPVECTOR is an alternate form of TwistExp()

    assert(isvector(xi), "This function only handles vectors");

    esp = TwistExp(HomogeneousToTwist(xi), theta); 
end