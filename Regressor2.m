function Y = Regressor2(DHtable, q, qp, v, vp, t, g0)
    n = size(DHtable,1);
    Y = Regressor(DHtable, q, qp, v, vp, t, g0, 1);
    
    for k = [2, n]
        Yk = Regressor(DHtable, q, qp, v, vp, t, g0, k);
        Y  = [Y, Yk];
    end
end