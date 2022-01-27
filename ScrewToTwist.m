function twist = ScrewToTwist(h, q, w)
% SCREWTOTWIST(h, q, w) builds a twist TWIST from its screw coordinates H 
%   (pitch), Q (point on axis), W (axis).

    twist = [-AxisToSkew(w) * q + h * w, w];
end