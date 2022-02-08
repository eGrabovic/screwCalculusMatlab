function out = rotXvectorized(x)
% ROTXVECTORIZED(x) returns a multidimensional array OUT which is composed
%   of rotation matrices. Each element of X is considered as an angle value
%   representing a rotation along the X-axis, from which every rotation
%   matrix in OUT is derived.

    C = cos(x);
    S = sin(x);
    
    [r, c] = size(C);
    
    Zer = zeros(r, c);
    Ones = ones(r, c);
    
    C    = reshape(C, 1, 1, r, c);
    S    = reshape(S, 1, 1, r, c);
    Zer  = reshape(Zer, 1, 1, r, c);
    Ones = reshape(Ones, 1, 1 ,r, c);
    
    out = [Ones, Zer,  Zer; 
           Zer,  C,   -S; 
           Zer,  S,    C];

end