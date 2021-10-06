function out = spatialTwistDerivative(spatialJac, qdot, qdotdot)
%
%
%
out = spatialJac*qdotdot + spatialJacDerivative(spatialJac, qdot)*qdot;

end