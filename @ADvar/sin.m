function adObj = sin(adObj)
% SIN(adObj) overloads the classic sin() function for the ADvar class.

    adObj.der = cos(adObj.val) .* adObj.der;
    adObj.val = sin(adObj.val);
    
end