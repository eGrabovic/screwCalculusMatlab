function res = ObjectJac(func, pars)
    res = [[eye(3),     zeros(3)]
           [zeros(3),   func(pars)]];
end