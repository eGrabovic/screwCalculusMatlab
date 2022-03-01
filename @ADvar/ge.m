function out = ge(adObj, comparison)
% GE(adObj, comparison) overloads the classic ge() function for the 
%   ADvar class.

    % TOOD: first if useless because this is a class method
    if ~isa(adObj,'ADvar')
        out = adObj >= comparison.val;
        return
    
    end
    
    if ~isa(comparison,'ADvar')
        out = adObj.val >= comparison;
        return
    
    end
    
    out = adObj.val >= comparison.val;

end