function s = squad(x0,x1,x2,x3,t)
% TODO: what is this?
    s = (1 - 2*t * (1-t)) * ((1-t) * x0 + t * x3) + ...
        2*t * (1-t) * ((1-t) * x1 + t * x2);
end
