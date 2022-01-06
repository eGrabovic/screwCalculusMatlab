function twist = RevoluteTwist(q,w)
% Gives Xi 6 vector given a point on axis and axis unit vector for a Revolute Joint *)
    c = [Cross(q,w), w];
    twist = reshape(c.',1,[]);    
end
