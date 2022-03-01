function adObj = atan2(adObj, X)
% ATAN2(adObj, X) overloads the classic atan2() function for the ADvar class.

    adObj.der = X.val * adObj.der / (X.val .^ 2 + adObj.val .^ 2) - ...
                adObj.val * X.der / (X.val .^ 2 + adObj.val .^ 2);
    adObj.val = atan2(adObj.val, X.val);

end