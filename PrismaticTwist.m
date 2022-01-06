function twist = PrismaticTwist(q, w)
% Gives Xi 6 vector given a point on axis and axis unit vector
% for a Prismatic Joint
    mat = [w, [0,0,0]];
    twist = reshape(mat.',1,[]);
end
