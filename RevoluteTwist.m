function twist = RevoluteTwist(q, w)
% REVOLUTETWIST(q,w) builds the 6-element vector TWIST corresponding to 
%   point Q on the axis with unit vector W for a prismatic joint

    c = [Cross(q,w), w];
    twist = reshape(c.',1,[]);    
end
