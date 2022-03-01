function adObj = sqrt(adObj)
% SQRT(adObj) overloads the classic sqrt() function for the ADvar class.

    adObj.der = 1 ./ (2 * sqrt(adObj.val)) .* adObj.der;
    adObj.val = sqrt(adObj.val);
    
end