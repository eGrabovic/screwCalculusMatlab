function ForceClose2DQFrames(qrows)
    lastrow = qrows(1);
    n = size(qros, 2);
    
    for i = [1, n]
        thisrow = ForceCloseQFrames(qrows(i));
        
        if thisrow(1) * lastrow(1) >= 0
            lastrow = thisrow;
        else
            lastrow = - thisrow;
        end
    end
end