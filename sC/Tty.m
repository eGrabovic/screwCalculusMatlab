function out = Tty(y)

out = eye(4, class(y));
out(2,4) = y;


end