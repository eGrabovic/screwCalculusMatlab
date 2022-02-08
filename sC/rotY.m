function R = rotY(alpha)
% ROTY(alpha) describes a rotation of the Euler angle ALPHA around
%   the Y-axis.
%   Note that the output matrix R is always with angles in radiants.
%   Allow for ADvar class argument ALPHA.
%
%   Examples TODO: expand with others
%       rotY(pi)
%   --> [-1   0   0
%         0   1   0
%         0   0  -1]
%
%   Input
%       alpha:      angle of rotation [rad or deg]
%   Output
%       R:          3x3 rotation matrix
%
    
    % Manage ADvar argument
    if isa(alpha, 'ADvar')
        T = rotY(alpha.val);
        Der = hat([0; 1; 0]) * T .* alpha.der;
        R = ADvar(T, Der);
        return
    end
    
    % Standard case
    
    R = [cos(alpha)  0  sin(alpha); 
         0           1  0; 
        -sin(alpha)  0  cos(alpha)];
end