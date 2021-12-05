function v = HomogeneousToVector(hom_v)
% HomogeneousToVector(...) extracts the cartesian coordinates v from
%   the homogeneous components HOM_V.
%
%   Example:
%       HomogeneousToVector([1.5, 2.5, 3.5, 1])
%   --> [1.5, 2.5, 3.5]
%
%   Input
%       hom_v:  numerical array, 4x1 or 1x4
%   Output
%       v:      numerical array, 3x1 or 1x3 (depending on input)
%
    assert(all(size(hom_v) == [1 4]) || all(size(hom_v) == [4 1]), ...
            "Homogeneous vectors must be of four elements.");

    if size(hom_v, 1) > size(hom_v,2) % column array
        v = hom_v(1:3,1);
    else % row array
        v = hom_v(1,1:3);
    end
end