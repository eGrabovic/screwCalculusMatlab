function R = diffRotZ2D(alpha)
% DIFFROTZ2D returns a R2 matrix R describing the differentiate of a Z-axis
%   rotation of angle ALPHA.

    % TODO: checks? alpha numeric?

    R = [[-sin(alpha),  -cos(alpha)]
         [cos(alpha),   -sin(alpha)]];

end