function pitch = TwistPitch(xi)
% TWISTPITCH(xi) gives the pitch PITCH of screw corresponding to a twist XI

% Find the pitch associated with a twist in R^3

    assert(isvector(xi), "xi is not a vector");
    assert(max(size(xi) == 6), "xi has wrong dimensions for a twist");
    
    c = mat2cell(xi, 1, [3 3]);
    v = c{1};
    w = c{2};
    
    pitch = (v * w) / (w * w);

end