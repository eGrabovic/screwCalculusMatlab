function mat = Hat(v)
    if (size(v,1) ~= 3)
        return
    end
    
    mat = [[ 0,       -v(3),        v(2)]
         [ v(3),     0,          -v(1)]
         [-v(2),    v(1),         0]];
end