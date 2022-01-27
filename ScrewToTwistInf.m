function twist = ScrewToTwistInf(q, w)
% SCREWTOTWISTINF(q, w) is an edge case of ScrewToTwist()

% Build a twist from a screw
    twist = [w, [0,0,0]];
end
