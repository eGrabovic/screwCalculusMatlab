function out = ad(twist, varargin)
%
%
%
w = twist(4:6);
v = twist(1:3);
w_hat = hat(w);
out = [w_hat, hat(v); zeros(3, 3), w_hat];

end