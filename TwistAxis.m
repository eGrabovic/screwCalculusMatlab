function [axis, angle] = TwistAxis(xi)
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