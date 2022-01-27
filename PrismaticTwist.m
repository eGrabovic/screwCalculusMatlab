function twist = PrismaticTwist(q, w)
% PRISMATICTWIST(q,w) builds the 6-elements vector TWIST corresponding 
%   to point Q on the axis with unit vector W for a prismatic joint

    mat = [w, [0,0,0]];
    twist = reshape(mat.',1,[]);
end
