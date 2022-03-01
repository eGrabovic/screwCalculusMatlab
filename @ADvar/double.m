function out = double(adObj)
% DOUBLE(adObj) overloads the classic double() function for the ADvar class.

    out = [adObj.val, adObj.der];
    
end