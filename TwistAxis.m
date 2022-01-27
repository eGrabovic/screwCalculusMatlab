function [axis, angle] = TwistAxis(xi)
% TWISTAXIS(xi) gives the axis AXIS of screw corresponding to a twist XI as
%   and the relative ANGLE
%   TODO: check if this is correct^

% Find the axis of a twist

    assert(isvector(xi), "xi is not a vector");

    c = mat2cell(xi, 1, [3 3]);
    v = c{1};
    w = c{2};

     if all(w = [0,0,0]) || all(w == [0; 0; 0])
        axis = 0;
        angle = v / sqrt(v.v);
     else 
         axis = (AxisToSkew(w) * v) / (w * w);
         angle = w / (w * w);
     end
end