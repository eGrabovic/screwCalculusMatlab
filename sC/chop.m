function var = chop(var, tol)
% CHOP(var, tol) replaces approximate real numbers in VAR that are close to
%   zero by the exact integer 0.
%   Uses a default tolerance of 10^-10.
%   Coded to mimic Wolfram Mathematica 'chop[]', even if in MATLAB it
%   already exists the 'round()' function.

    if exist('tol', 'var') == 0
        tol = 1e-10;
    end

    var = round(var./tol).*tol;
end