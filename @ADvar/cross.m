function adObj = cross(adObj, multiplier)
% CROSS(adObj, multiplier) overloads the classic cross() function for the 
%   ADvar class.
            
   adObj.der = cross(adObj.der, multiplier.val) + cross(adObj.val, multiplier.der);
   adObj.val = cross(adObj.val, multiplier.val);
    
end