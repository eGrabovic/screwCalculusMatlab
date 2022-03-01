function out = minus(adObj, addend)
% MINUS(adObj, addend) overloads the classic minus() function for the 
%   ADvar class.

    % TODO: first if useless because this is a class method
    if ~isa(adObj,'ADvar')
        addend.der = -addend.der;
        addend.val = adObj - addend.val;
        out = addend;

    elseif ~isa(addend,'ADvar')
        adObj.der = adObj.der;
        adObj.val = adObj.val - addend;
        out = adObj;

    else
        adObj.der = adObj.der - addend.der;
        adObj.val = adObj.val - addend.val;
        out = adObj;
        
    end
end