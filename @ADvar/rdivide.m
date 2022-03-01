function out = rdivide(adObj, divisor)
% RDIVIDE(adObj, divisor) overloads the classic rdivide() function for the 
%   ADvar class.
         
    % TODO: first if useless because this is a class method
    if isa(adObj,'ADvar') == false
        divisor.der = -adObj .* divisor.der ./ divisor .^ 2;
        divisor.val = adObj ./ divisor.val;
        out = divisor;
        
    elseif isa(divisor,'ADvar') == false
        adObj.der = adObj.der ./ divisor;
        adObj.val = adObj.val ./ divisor;
        out = adObj;
        
    else
        adObj.der = (adObj.der .* divisor.val - adObj.val .* divisor.der) ./ ...
                    (divisor.val .^ 2);
        adObj.val = adObj.val ./ divisor.val;
        out = adObj;
        
    end
end