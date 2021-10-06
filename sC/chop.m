function var = chop(var, tol)
var = round(var./tol).*tol;
end