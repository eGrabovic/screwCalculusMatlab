function mag = TwistMagnitude(xi)
% TWISTMAGNITUDE(xi) gives magnitude MAG of screw corresponding to a twist
%   XI.
  
    assert(isvector(xi), "xi must be a vector");

    [v, w] = Partition(xi, 3);

    if all(w == [0,0,0]) || all(w == [0; 0; 0]) 
      mag = sqrt(v * v);
    else
        mag = sqrt(w * w);
    end
end