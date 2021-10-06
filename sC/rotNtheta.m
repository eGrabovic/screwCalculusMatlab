function out = rotNtheta(n, theta)

out = n*n.' + (eye(3) - n*n.').*cos(theta) + hat(n).*sin(theta);

end