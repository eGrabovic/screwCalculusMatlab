function res = SkewExp(S, theta)
% Matrix exponential for a skew symmetric matrix or a vector

    assert(isSkew(S), "Provided matrix must be a skew matrix");
    assert(isnumeric(theta), "Provided angle must be a rumeric value");

    % Use Rodrigues' formula
    res = eye(3) + sin(theta) * S + (1 - cos(theta)) * S * S;
end