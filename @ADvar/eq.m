function out = eq(adObj, comparison)
% EQ(adObj, comparison) overloads the classic eq() function for the 
%   ADvar class.

    % TOOD: first if useless because this is a class method
    if ~isa(adObj,'ADvar')
        out = adObj == comparison.val;
        return
    
    end
    
    if ~isa(comparison,'ADvar')
        out = adObj.val == comparison;
        return
    
    end
    
    out = adObj.val == comparison.val;

end