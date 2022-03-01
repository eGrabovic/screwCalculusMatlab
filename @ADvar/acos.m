function adObj = acos(adObj)
% ACOS(adObj) overloads the classic acos() function for the ADvar class.

    adObj.der = -1 ./ sqrt(1 - adObj.val * adObj.val) * adObj.der;
    adObj.val = acos(adObj.val);
    
end