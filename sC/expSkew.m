function X = expSkew(axis, th)
% EXPSKEW(axis, th) computes the exponenxial X such that 
%   X = e^(axis * th)
%   with 
%   - AXIS: axis versor
%   - TH:  angle
%
%   This uses Rodriguez' formula if the supplied axis is generic, and
%   copies any elementary rotation if the axis coincides with one of the
%   main axis (to improve performance).

    if all(isnumeric(axis)) == true
        if abs(axis.'*axis) -1 >= 1e-8
            error('axis requires to be a versor');
        end
        if abs(axis(1)) <= 1e-8 && abs(axis(2)) <= 1e-8
            
            if axis(3) > 0
                
                X = rotZ(th);
                
            elseif axis(3) < 0
                
                X = rotZ(-th);
            end
            
        elseif abs(axis(2)) <= 1e-8 && abs(axis(3)) <= 1e-8
            
            if axis(1) > 0
                
                X = rotX(th);
                
            elseif axis(1) < 0
                
                X = rotX(-th);
            end
            
        elseif abs(axis(1)) <= 1e-8 && abs(axis(3)) <= 1e-8
            
            if axis(2) > 0
                
                X = rotY(th);
                
            elseif axis(2) < 0
                
                X = rotY(-th);
            end
        else
            axisHat = hat(axis);
            X = eye(3) + axisHat*sin(th) + axisHat*axisHat*(1 - cos(th));
            %rodriguez formula
        end
        
    else
        if abs(double(axis(1))) <= 1e-8 && abs(double(axis(2))) <= 1e-8
            
            if double(axis(3)) > 0
                X = rotZ(th);
            elseif double(axis(3)) < 0
                X = rotZ(-th);
            end
        elseif abs(double(axis(1))) <= 1e-8 && abs(double(axis(3))) <= 1e-8
            
            if double(axis(2)) > 0
                X = rotY(th);
            elseif double(axis(2)) < 0
                X = rotY(-th);
            end
        elseif abs(double(axis(2))) <= 1e-8 && abs(double(axis(3))) <= 1e-8
            
            if double(axis(1)) > 0
                X = rotX(th);
            elseif double(axis(1)) < 0
                X = rotX(-th);
            end
        else
            X = eye(3) + hat(axis)*sin(th) + hat(axis)*hat(axis)*(1 - cos(th));
            % Rodriguez formula
        end
    end
end