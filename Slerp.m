function r = Slerp(p0, p1, t)
% SLERP(p0, p1, t) returns the Spherical Linear intERPolation function R
%   bewteen the unit quaternions P0 and P1.
%   Slerp(p0, p1, t=0) = p0 
%   and 
%   Slerp(p0, p1, t=1) = p1

    costh = round(p0 * p1);
    
    if costh > 0.0
        costh = round(costh - 1) + 1;
    end
    if costh < 0.0
        costh = round(costh + 1) - 1;
    end

    th = acos(costh);
    sinth = sin(th);
    
    if sinth == 0
        r = (1-t)*p0 + t*p1;
    else
        r = (sin(th*(1-t)) / sinth) * p0 + (sin(th*t) / sinth) * p1;
    end
end