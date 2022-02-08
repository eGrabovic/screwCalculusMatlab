function vet = SkewToAxis(mat)
% SKEWTOAXIS(mat) generates an axis VET from a skew symmetric matrix MAT

    % TODO: implement skewQ and use it here
    assert(all(size(mat) == [3, 3]), ...
        "Argument has wrong dimensions!"); % and not skewQ

    vet = [mat(3,2), mat(1,3), mat(2,1)];
end