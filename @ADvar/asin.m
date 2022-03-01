function adObj = asin(adObj)
% ASIN(adObj) overloads the classic asin() function for the ADvar class.

    adObj.der = 1 ./ sqrt(1 - adObj.val .* adObj.val) .* adObj.der;
    adObj.val = asin(adObj.val);
    
end