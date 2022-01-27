function out = ad(twist, varargin)
% AD(TWIST, VARARGIN) gives the adjoint matrix OUT corresponding to the
%   twist TWIST.
%   TODO: varargin not needed
  
    w = twist(4:6);
    v = twist(1:3);

    w_hat = hat(w);

    out = [w_hat, hat(v); 
           zeros(3, 3), w_hat];

end