function res = SkewExpVect(v, theta)
    assert(isvector(v), "v must be a vector");
    res = SkewExp(AxisToSkew(v),theta);
end
