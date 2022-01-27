function G = GlobalGraspMatrix(EPlist)
% GLOBALGRASPMATRIX(EPlist) returns the global grasp matrix G for the 
%   changes of reference point expressed by the vector list EPLIST.
%   Beware that the EPs are the components of the vector EP in a fixed 
%   reference frame.

    G = [];
    
    for i = [1, size(EPlist,2)]
        G = [G, GraspMatrix(EPlist(i))];
    end
end