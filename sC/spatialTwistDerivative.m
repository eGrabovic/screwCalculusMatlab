function out = spatialTwistDerivative(spatialJac, qdot, qdotdot)
% SPATIALTWISTDERIVATIVE(jac, qdot) computes the derivate OUT of the 
%   tewist of a robot described through the spatial Jacobian JAC w.r.t. 
%   the joint variables derivate QDOT and QDOTDOT.

    out = spatialJac * qdotdot + ...
          spatialJacDerivative(spatialJac, qdot) * qdot;

end