function [axis, theta] = RotationParam(Rotmat)
    [nr, nc] = size(Rotmat);

    % Check to make sure that our input makes sense
    assert(ismatrix(Rotmat) && (nr == nc) && (nr == 3), ...
        "Wrong matrix dimensions");     
     
    t = (trace(Rotmat) - 1) / 2;

    if (t < -1)
        t = -1;
    elseif (t > 1)
        t = 1;
    end

    theta = acos(t);

    if (theta ~= 0)
       axis = RotationAxis(Rotmat, theta);
    else
       axis = Table(0, {nc});
    end

end