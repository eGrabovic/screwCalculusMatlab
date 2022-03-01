function out = plus(adObj, addend)
% PLUS(adObj, addend) overloads the classic plus() function for the 
%   ADvar class.

    % TODO: first if useless because this is a class method
    if isa(adObj,'ADvar') == false
        v.der = v.der;
        v.val = adObj + v.val;
        out = v;

    elseif isa(v,'ADvar') == false
        adObj.der = adObj.der;
        adObj.val = v + adObj.val;
        out = adObj;

    else
        adObj.der = adObj.der + v.der;
        adObj.val = adObj.val + v.val;
        out = adObj;
        
    end
end