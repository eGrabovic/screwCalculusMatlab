function adObj = atan(adObj)
% ATAN(adObj) overloads the classic atan() function for the ADvar class.

    adObj.der = adObj.der ./ (adObj.val .^ 2 + 1);
    adObj.val = atan(adObj.val);
    
end