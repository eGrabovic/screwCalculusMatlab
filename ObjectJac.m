function JAC = ObjectJac(func, pars)
% OBJECTHAC(function, pars) return the 6x6 object jacobian JAC 
%   corresponding to the parametrization function FUNC applied to 
%   parameters PARS

    JAC = [[eye(3),     zeros(3)]
           [zeros(3),   func(pars)]];
end