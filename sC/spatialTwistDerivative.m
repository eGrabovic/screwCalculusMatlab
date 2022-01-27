function out = spatialTwistDerivative(spatialJac, qdot, qdotdot)
% TODO

out = spatialJac*qdotdot + spatialJacDerivative(spatialJac, qdot)*qdot;

end