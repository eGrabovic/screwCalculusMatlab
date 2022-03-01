function adObj = cos(adObj)
% COS(adObj) overloads the classic cos() function for the ADvar class.

    adObj.der = -sin(adObj.val) .* adObj.der;
    adObj.val = cos(adObj.val);
    
end