function ForceClose3DQFrames(qplanes)
% TODO: unspecified
    lastplane = qplanes(1);
    n = size(qplanes, 2);
    
    for i = [1, n]
        thisplane = ForceClose2DQFrames(qplanes(i));
        
        if thisplane(1) * lastplane(1) >= 0
            lastplane = thisplane;
        else
            lastplane = -thisplane;
        end
    end
end