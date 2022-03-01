function adObj = log(adObj)
% LOG(adObj) overloads the classic log() function for the ADvar class.

    adObj.der = adObj.der ./ adObj.val;
    adObj.val = log(adObj.val);
    
end