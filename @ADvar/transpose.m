function adObj = transpose(adObj)
% TRANSPOSE(adObj) computes the transpose of each property of the class 
%   object ADOBJ.

    adObj.der = adObj.der.';
    adObj.val = adObj.val.';
    
end