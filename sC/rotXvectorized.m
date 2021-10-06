function out = rotXvectorized(x)

C = cos(x);
S = sin(x);

[r, c] = size(C);

Zer = zeros(r, c);
Ones = ones(r, c);

C = reshape(C, 1, 1, r, c);
S = reshape(S, 1, 1, r, c);
Zer = reshape(Zer, 1, 1, r, c);
Ones = reshape(Ones, 1, 1 ,r, c);

out = [Ones, Zer, Zer; Zer, C, -S; Zer, S, C];



end