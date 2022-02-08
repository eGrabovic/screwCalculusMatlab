function R = rotX(alpha)
% ROTX(alpha) describes a rotation of the Euler angle ALPHA around
%   the X-axis.
%   Note that the output matrix R is always with angles in radiants.
%   Allow for ADvar class argument ALPHA.
%
%   Examples TODO: expand with others
%       rotX(pi)
%   --> [1   0   0
%        0  -1   0
%        0   0  -1]
%
%   Input
%       alpha:      angle of rotation [rad or deg]
%   Output
%       R:          3x3 rotation matrix
%

    % Manage ADvar argument
    if isa(alpha, 'ADvar')
        T = rotX(alpha.val);
        Der = hat([1; 0; 0]) * T .* alpha.der;
        R = ADvar(T, Der);
        return
    end

    % Standard case
    
    R = [1  0           0;
         0  cos(alpha)  -sin(alpha);
         0  sin(alpha)  cos(alpha)];

end