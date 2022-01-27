function len = QLength(curveOnSphere)
% QLENGTH(curveOnSphere) returns ...
%   TODO: what is this exactly? distance between what?

% Distances in quaternion arcs
    n = size(curveOnSphere, 2) - 1;
    len = [];

    for i = [1, n]
        len = len + acos(curveOnSphere(i) . curveOnSphere(i+1));
    end
end