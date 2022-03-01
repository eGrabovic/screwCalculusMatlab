function adObj = tan(adObj)
% TAN(adObj) overloads the classic tan() function for the ADvar class.

    adObj.der = 2 ./ (cos(2 .* adObj.val) + 1);
    adObj.val = tan(adObj.val);
    
end