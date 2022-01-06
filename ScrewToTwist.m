function twist = ScrewToTwist(h, q, w)
    twist = [-AxisToSkew(w) * q + h * w, w];
end