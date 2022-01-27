function axis = RotationAxis(Rotmat, theta)
% ROTATIONAXIS(Rotmat, theta) returns the rotation axis AXIS of ROTMAT 
%   in SO(3)  

    [nr, nc] = size(Rotmat);

    % Check to make sure that our input makes sense
    assert(ismatrix(Rotmat) && (nr == nc) && (nr == 3), ...
        "Wrong matrix dimensions");

    asert((theta >= 0) && (theta <= pi), ...
           "Theta is not in a valid interval");
 
    if (theta == pi)
      nullspace = null(Rotmat - eye(3));
      axis = nullspace(1);
      axis = axis / norm(axis);
    else
      axis = [Rotmat(3, 2) - Rotmat(2, 3)
              Rotmat(1, 3) - Rotmat(3, 1)
              Rotmat(2, 1) - Rotmat(1, 2)] / (2 * sin(theta));
    end
end