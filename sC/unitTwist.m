function uT = unitTwist(jointType, axis, q)
% UNITTWIST(jointType, axis, q) creates a unit twist UT from 
%   - the type of joint JOINTTYPE ('P' for a prismatic, 'R' for a rotoidal)
%   - the joint axis AXIS
%   - the distance Q from the origin in case of a 'R' joint

    % Prismatic joint
    if strcmpi(jointType,'P') == true
        uT = [axis; 0; 0; 0];
    end
    
    % Rotoidal joint
    if strcmpi(jointType,'R') == true
        uT = [-cross(axis,q); axis];
    end

end