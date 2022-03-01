function out = le(adObj, comparison)
% LE(adObj, comparison) overloads the classic le() function for the 
%   ADvar class.
     
    % TODO: first if useless because this is a class method 
    if ~isa(adObj,'ADvar')
        out = adObj <= comparison.val;
        return
    
    end
    
    if ~isa(comparison,'ADvar')
        out = adObj.val <= comparison;
        return
    
    end
    
    out = adObj.val <= comparison.val;
     
 end