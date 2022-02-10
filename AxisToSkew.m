function skew = AxisToSkew(omega)
% AXISTOSKEW(omega) generates a skew symmetric matrix SKEW from the axis 
%   described through the array OMEGA.

    % Check to make sure the dimensions are okay
    if(length(omega)~=3)
        disp("Argument has wrong dimensions!");
        return
    end

    % Return the appropriate matrix
    skew = [[0,         -omega(3),  omega(2)  ]
            [ omega(3), 0,          -omega(1) ]
            [-omega(2), omega(1),   0         ]];
end