function prod = QuatProd(q, p)
% QUATPROD(q, p) returns the composition product PROD between two 
%   quaternions Q and P

% TODO: check arguments
prod = [ q(1)*p(1) - q(2)*p(2) - q(3)*p(3) - q(4)*p(4), ...
         q(1)*p(2) + q(2)*p(1) + q(3)*p(4) - q(4)*p(3), ...
         q(1)*p(3) + q(3)*p(1) + q(4)*p(2) - q(2)*p(4), ...
         q(1)*p(4) + q(4)*p(1) + q(2)*p(3) - q(3)*p(2)];
end