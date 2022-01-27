function G = GraspMatrix(EP)
% GRASPMATRIX(EP) returns the grasp matrix G for a change of reference 
%   point expressed by the vector EP.
%   Beware that the components of EP are the components of the vector EP 
%   in a fixed reference frame, not in the local E.E. frame.

    G = [[eye(3),   zeros(3,3)];
        [Hat(EP),   eye(3)]];
end