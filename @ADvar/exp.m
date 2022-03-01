function adObj = exp(adObj)
% EXP(adObj) overloads the classic exp() function for the ADvar class.

    adObj.der = adObj.der .* exp(adObj.val);
    adObj.val = exp(adObj.val);
    
end