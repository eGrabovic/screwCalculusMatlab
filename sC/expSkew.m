function Xexp = expSkew(axis, th)
% EXPSKEW(axis, th) computes the exponenxial XEXP such that 
%   Xexp = e^(axis * th)
%   with 
%   - AXIS: axis versor
%   - TH:   angle
%
%   This uses Rodriguez' formula if the supplied axis is generic, and
%   copies any elementary rotation if the axis coincides with one of the
%   main axis (to improve performance).

% TODO: how to add formulas?

    assert(isnumeric(theta), "Provided angle must be a rumeric value");

    threshold = 1e-8;

    if all(isnumeric(axis)) == true

        % Check if axis is a unit vector within an acceptable threshold
        if abs(axis.'*axis) -1 >= threshold % Square of norm
            error('axis must be a versor');
        end

        abs_x = abs(axis(1));
        abs_y = abs(axis(2));
        abs_z = abs(axis(3));
    
        % Check if the Z-axis is specified, within a threshold
        if checkZeroThreshold([abs_x, abs_y], threshold)
            % Depending on the angle sign
            Xexp = rotZ((2*sign(axis(3))-1) * th);
            
        % Check if the X-axis is specified, within a threshold
        elseif checkZeroThreshold([abs_y, abs_z], threshold)
            % Depending on the angle sign
            Xexp = rotX((2*sign(axis(1))-1) * th);
            
        % Check if the Y-axis is specified, within a threshold
        elseif checkZeroThreshold([abs_x, abs_z], threshold)
            % Depending on the angle sign
            Xexp = rotY((2*sign(axis(2))-1) * th);

        else % no particular axis has been provided
            % Use Rodriguez' formula
            axisHat = hat(axis);
            Xexp = eye(3) + axisHat * sin(th) + ...
                axisHat * axisHat * (1 - cos(th));
        end
        
    else % TODO: which types to consider?

        % TODO: check if data can be converted to double

        abs_x = abs(double(axis(1)));
        abs_y = abs(double(axis(2)));
        abs_z = abs(double(axis(3)));

        if checkZeroThreshold([abs_x, abs_y], threshold)
            % Depending on the angle sign
            Xexp = rotZ((2*sign(double(axis(3))-1)) * th);

        elseif checkZeroThreshold([abs_x, abs_z], threshold)
            % Depending on the angle sign
            Xexp = rotY((2*sign(double(axis(2))-1)) * th);

        elseif checkZeroThreshold([abs_y, abs_z], threshold)
            % Depending on the angle sign
            Xexp = rotX((2*sign(double(axis(1))-1)) * th);

        else % no particular axis has been provided
            % Rodriguez formula
            Xexp = eye(3) + hat(axis) * sin(th) + ...
                hat(axis) * hat(axis) * (1 - cos(th));
        end
    end
end