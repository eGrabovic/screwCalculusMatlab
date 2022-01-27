function lastfrm = ForceCloseQFrames(qfrms)
% TODO: unspecified
    lastfrm = qfrms(1);
    n = size(qfrms, 2);
    
    for i = [1, n]
        thisfrm = qfrms(i);
        
        if thisfrm * lastfrm >= 0
            lastfrm = thisfrm;
        else
            lastfrm = - thisfrm;
        end
    end
end