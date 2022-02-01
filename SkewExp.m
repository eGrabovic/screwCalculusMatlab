function res = SkewExp(S, theta)
% SKEWEXP(S,theta) gives the matrix exponential RES of an axis S.
%   Default value of THETA is 1.

% Matrix exponential for a skew symmetric matrix or a vector

% TODO: this is an alternative to expSkew that directly uses a Skew matrix;
%   change name?

    assert(isSkew(S), "Provided matrix must be a skew matrix");
    assert(isnumeric(theta), "Provided angle must be a rumeric value");

    % Use Rodrigues' formula
    res = eye(3) + sin(theta) * S + (1 - cos(theta)) * S * S;
end