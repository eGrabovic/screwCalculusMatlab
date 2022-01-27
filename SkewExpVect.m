function res = SkewExpVect(v, theta)
% SKEWEXPVECT(v, theta) is an edge case of SkewExp()

    assert(isvector(v), "v must be a vector");
    res = SkewExp(AxisToSkew(v),theta);
end
