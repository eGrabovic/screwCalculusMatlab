function r = Slerp(p0,p1,t)
    costh = round(p0 . p1);
    
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