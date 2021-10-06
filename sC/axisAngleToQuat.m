function out = axisAngleToQuat(theta, n)
out = zeros(4,1);
out(1) = cos(theta/2);
out(2) = sin(theta/2).*n(1);
out(3) = sin(theta/2).*n(2);
out(4) = sin(theta/2).*n(3);
end