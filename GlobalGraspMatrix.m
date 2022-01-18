function G = GlobalGraspMatrix(EPlist)
    G = [];
    
    for i = [1, size(EPlist,2)]
        G = [G, GraspMatrix(EPlist(i))];
    end
end