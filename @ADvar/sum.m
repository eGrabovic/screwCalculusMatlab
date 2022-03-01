function adObj = sum(adObj, dim)
% SUM(adObj, dim) returns the ADvar object whose properties are the 
%   sum of the property of the ADvar object passed as argument, along the 
%   dimension DIM.

% TODO: why the dimension? aren't the properties integers?
            
    adObj.val = sum(adObj.val, dim);
    adObj.der = sum(adObj.der, dim);

end