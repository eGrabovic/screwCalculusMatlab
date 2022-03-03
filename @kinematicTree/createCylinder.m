function [x, y, z] = createCylinder(r, len)
% TODO

% TODO: why as method?
    theta = linspace(0, 2 * pi, 15);
    x = r .* cos(theta);
    y = r .* sin(theta);
    z = len / 2 .* ones(1, 15);
    
end