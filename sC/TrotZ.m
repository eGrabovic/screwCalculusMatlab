function out = TrotZ(alpha)
% TROTZ(alpha) computes the homogeneous matrix OUT describing a rotation of
%   angle ALPHA around the z-axis.
%   Allows for ADvar arguments.

    % Handles ADvar argument
    if isa(alpha, 'ADvar')
        T = rotZ(alpha.val);
        Der = hat([0; 0; 1]) *T .* alpha.der;
        Der = [[Der, [0; 0; 0]]; % TODO: ok 0s as a column?
               [0,    0, 0, 0]];
        out = ADvar(TrotZ(alpha.val), Der);

        return
    end
    
    % Standard case
    out = [[cos(alpha), -sin(alpha),    0,  0];
           [sin(alpha),  cos(alpha),    0,  0];
           [0,           0,             1,  1]];

end