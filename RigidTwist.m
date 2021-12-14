function [twist, theta] = RigidTwist(g)
% Find the twist which generates a rigid motion
    
    assert(ismatrix(g), "g is not a matrix");

    % Make sure the dimensions are okay
    % ! Missing ! TODO

    % Extract the appropriate pieces of the homogeneous transformation
    Rot = RigidOrientation(g);
    pos = RigidPosition(g);

    % Now find the axis from the rotation
    [w, theta] = RotationParam(Rot);
    hatw = Hat(w);
     
    % Split into cases depending on whether theta is zero *)
    if (theta == 0)
        theta = norm(pos);
        if (theta == 0)
            twist = [0,0,0,0,0,0];
            theta = 0;
        else
            v = pos/theta;
        end
    else
        % Solve a linear equation to figure out what v is
        Ainv = eye(3)/theta - (1/2) * hatw + ...
                        (1/theta - (1/2) * cot(theta/2)) * hatw.^2; 
        v = Ainv.p;
    end

    twist = [v, w];
    % theta already assigned
	
end