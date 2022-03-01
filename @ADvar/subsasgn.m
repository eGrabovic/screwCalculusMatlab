function adObj = subsasgn(adObj, s, val)
% See the documentation for subsasgn()

% TODO ^

    switch s(1).type
        case '.'
            adObj = builtin('subsasgn', adObj, s, val);
        case '()'
            if isa(val, 'ADvar')
                adObj.der = subsasgn(adObj.der, s, val.der);
                adObj.val = subsasgn(adObj.val, s, val.val);
            else
                adObj.der = subsasgn(adObj.der, s, 0);
                adObj.val = subsasgn(adObj.val, s, val);
            end
    end
end