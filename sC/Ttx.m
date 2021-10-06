function out = Ttx(x)

out = eye(4, class(x));
out(1,4) = x;

end