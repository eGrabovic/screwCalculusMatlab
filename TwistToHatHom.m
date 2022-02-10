function X = TwistToHatHom(vec)
% HAT(vec) transforms a R3 column twist vector VEC in a homogeneous 'hat' 
%   R4x4 matrix X.
   
    assert(isvector(vec), ...
           "Provided argument must be an array");

    X = [[hat([vec(4); vec(5); vec(6)]),  vec(1:3)]
         [0,   0,   0,                    0]];

end