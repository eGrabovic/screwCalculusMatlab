function quat = MakeQuat(n, angle)
% MAKEQUAT(axis, angle) creates the unit quaternion QUAT given the unit 
%   vector N of the axis and the ANGLE

% WARNING: these are the reverse of conventional argument order;
% change order once all is consistent to (angle, nhat). TODO?

    c = cos(angle/2);
    s = sin(angle/2);
    
    quat = [c, n(1)*s, n(2)*s, n(3)*s];
end