function out = TrotX(x)
% TROTX(x) computes the homogeneous matrix OUT describing a rotation of
%   angle X around the x-axis.
%   Allows for ADvar arguments.

    % Handles ADvar argument
    if isa(x, 'ADvar')
        T = rotX(x.val);
        Der = hat([1; 0; 0]) * T .* x.der;
        Der = [[Der, [0; 0; 0]]; % TODO: ok 0s as a column?
               [0,    0, 0, 0]];
        out = ADvar(TrotX(x.val), Der);

        return
    end
    
    % Standard case
    out = [[1,      0,       0,          0];
           [0,   cos(x),     -sin(x),    0];
           [0,   sin(x),     cos(x),     0];
           [0,   0,          0,          1]];

end