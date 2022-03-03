function g = jointFun(x, h)
% TODO

% TODO: why as a method?
    switch h
        case 0
            g = [[cos(x), -sin(x), 0,    0];
                 [sin(x),  cos(x), 0,    0];
                 [0,       0,      1,    0];
                 [0,       0,      0,    1]];
        case Inf
            g = [[eye(3),   [0; 0; x]]; 
                 [0, 0, 0,  1]];
        otherwise
            g = [[cos(x), -sin(x),  0,     0];
                 [sin(x),  cos(x),  0,     0];
                 [0,       0,       1,     x .* h];
                 [0,       0,       0,     1]];
    end
    
end