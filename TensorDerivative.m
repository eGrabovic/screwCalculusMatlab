function res = TensorDerivative(J, q)
% TENSORDERIVATIVE(J, q) computes the derivative RES of the (nxn) matrix J 
%   w.r.t. vector Q

    DJ  = diff(J, q);
    res =(1/2) * (DJ + DJ');
end