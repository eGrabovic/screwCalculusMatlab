function out = TrotY(x)
% TROTY(x) computes the homogeneous matrix OUT describing a rotation of
%   angle X around the y-axis.
%   Allows for ADvar arguments.

    % Handles ADvar argument
    if isa(x, 'ADvar')
        T = rotY(x.val);
        Der = hat([0; 1; 0]) *T .* x.der;
        Der = [[Der, [0; 0; 0]]; % TODO: ok 0s as a column?
               [0,    0, 0, 0]];
        out = ADvar(TrotY(x.val), Der);

        return
    end
    
    % Standard case
    out = [[cos(x),     0,  sin(x), 0];
           [0,          1,  0,      0];
           [-sin(x),    0,  cos(x), 0];
           [0,          0,  0,      1]];

end