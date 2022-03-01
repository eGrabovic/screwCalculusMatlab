function adObj = subsref(adObj, s)
% See the documentation for subsref()

% TODO ^

    switch s(1).type
        case '.'
            adObj = builtin('subsref', adObj, s);
        case '()'
            adObj.der = subsref(adObj.der, s);
            adObj.val = subsref(adObj.val, s);
    end

end