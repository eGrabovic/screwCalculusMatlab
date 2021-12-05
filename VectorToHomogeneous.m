function hom_v = VectorToHomogeneous(v)
% Convert a point into homogeneous coordinates

    % Check to make sure the dimensions of the args make sense
    assert(isvector(v), "The provided argument is not an array");

    % 1x1 scalar 'arrays' are not a vector
    assert(size(v,1) == size(v,2), ...
        "A single scalar has been provided; please provide an array");

    assert(size(v,1) <= 3 && size(v,2) <= 3, ...
        "Provided array is too long to represent a vector");

    % Now put everything together into a homogeneous array
    if size(a,1) > size(a,2) % column vector
        hom_v = [v; 0];
    else % row vector
        hom_v = [v, 0];
    end
end