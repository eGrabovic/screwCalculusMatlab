function var = chop(var, tol)
% TODO: Mathematica chop[]?
var = round(var./tol).*tol;
end