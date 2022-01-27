function hom_p = PointToHomogeneous(p)
% POINTTOHOMOGENEOUS(q) gives the homogeneous representation HOM_P of a 
%   point P

    % Check to make sure the dimensions of the args make sense
    assert(isvector(p), "The provided argument is not an array");

    % 1x1 scalar 'arrays' are not a point
    assert(size(p,1) == size(p,2), ...
        "A single scalar has been provided; please provide an array");

    assert(size(p,1) <= 3 && size(p,2) <= 3, ...
        "Provided array is too long to represent a point");

    % Now put everything together into a homogeneous point
    if size(a,1) > size(a,2) % column vector
        hom_p = [p; 1];
    else % row vector
        hom_p = [p, 1];
    end
end