function r = Lerp(p0, p1, t)
% LERP(p0, p1, t) returns the Linear intERPolation function between the 
%   unit quaternions P0 (initial) and P1 (final) for the time T.
%   Lerp(p0, p1, t=0) = p0
%   and
%   Lerp(p0, p1, t=1) = p1
% TODO: how to add formula?

    r = (1-t) * p0 + t * p1;
end