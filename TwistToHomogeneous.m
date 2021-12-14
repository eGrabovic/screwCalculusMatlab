function hom = TwistToHomogeneous(xi)
% Convert a twist to homogeneous coordinates

    assert(isvector(xi), "xi is not a vector");

    w = xitow(xi);
    v = xitov(xi);
    
    % Make sure that we got a real twist
    assert(all(size(w) ~= 0) && all(size(v) ~= 0), ...
        "Provided xi is not a real twist");

    % Now put everything together into a homogeneous transformation
    hom = [[AxisToSkew(w),  v']
           [zeros(1,3),     0 ];
end
