function out = times(adObj, multiplier)
% TIMS(adObj, multiplier) overloads the classic times() function for the 
%   ADvar class.

    % TODO: first if useless because this is a class method
    if ~isa(adObj,'ADvar')
        multiplier.der = adObj .* multiplier.der;
        multiplier.val = adObj .* multiplier.val;
        out = multiplier;
        
    elseif ~isa(multiplier,'ADvar')
        adObj.der = adObj.der .* multiplier;
        adObj.val = adObj.val .* multiplier;
        out = adObj;
        
    else
        adObj.der = adObj.der .* multiplier.val + adObj.val .* multiplier.der;
        adObj.val = adObj.val .* multiplier.val;
        out = adObj;
        
    end
end