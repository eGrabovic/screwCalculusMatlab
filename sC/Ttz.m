function out = Ttz(z)

out = eye(4, class(z));
out(3,4) = z;

end