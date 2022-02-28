function E = expTw(screwObj, unitTwist,th)
% EXPTW(unitTwist, th) computes the exponential E of a twist UNITTWIST w.r.t.
%   the joint variable TH.

    threshold = 1e-6;

    if all(isnumeric(unitTwist)) == true
        % Prismatic joint
        if checkZeroThreshold(unitTwist(4:6), threshold)
            axis = [unitTwist(1); unitTwist(2); unitTwist(3)]; % TODO: xitov?
            R = eye(3);
            d = axis*th;
            
        else % Revolute joint
            axis = [unitTwist(4); unitTwist(5); unitTwist(6)]; % TODO: xitow?
            v = [unitTwist(1); unitTwist(2); unitTwist(3)]; % TODO: xitov?

            q = -cross(v, axis);
            R = screwObj.expSkew(axis, th);
            d = (eye(3) - R) * q;
            
        end

        E = screwObj.hom_mat(R, d);

    end

end