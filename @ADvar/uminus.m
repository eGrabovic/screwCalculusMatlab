function adObj = uminus(adObj)
% UMINUS(adObj) overloads the classic uminus() function for the  ADvar class.
            
    adObj.val = -adObj.val;
    adObj.der = -adObj.der;
    
end