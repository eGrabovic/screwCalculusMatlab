function E = expTw(screwObj, unitTwist,th)
% EXPTW(unitTwist, th) computes the exponential E of a twist UNITTWIST w.r.t.
%   the joint variable TH.
%   Allows for ADvar class TH argument.
    
    % Handle th as an ADvar instance
    if isa(th, 'ADvar')
       E = expTw(unitTwist, th.val);
       E = ADvar(E, hat(unitTwist)*E.*th.der);
       return
    end

    threshold = 1e-6;

    % Standard case
    if all(isnumeric(unitTwist)) == true
        % Prismatic joint
        if checkZeroThreshold(unitTwist(4:6), threshold)
            % Extract translational components
            axis = xitov(unitTwist);
            
            % Make sure that we got a real twist
            assert(all(size(axis) ~= 0), ...
                   "The provided vector is not a real twist");

            R = eye(3);
            d = axis * th;
            
        else % Revolute joint
            % Extract angular and translational components
            axis = xitow(unitTwist);
            v = xitov(unitTwist);

            % Make sure that we got a real twist
            assert(all(size(axis) ~= 0) && all(size(v) ~= 0), ...
                   "The provided vector is not a real twist");

            q = -cross(v, axis);
            R = screwObj.expSkew(axis, th);
            d = (eye(3) - R) * q; % TODO: + (axis*(axis'*v)*th?
            
        end

        E = screwObj.hom_mat(R, d);

    end

end