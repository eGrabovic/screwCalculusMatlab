function adObj = abs(adObj)
% ABS(adObj) overloads the classic abs() function for the ADvar class.
    
    adObj.der = adObj.der * sign(adObj.val);
    adObj.val = abs(adObj.val);
    
end