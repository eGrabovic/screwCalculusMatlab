function out = power(adObj, pow)
% POWER(adObj, pow) overloads the classic power() function for the ADvar class.

    % TODO: first if useless because this is a class method
    if isa(adObj,'ADvar') == false
        pow.der = adObj .^ pow.val .* log(adObj) .* pow.der;
        pow.val = adObj .^ pow.val;
        out = pow;
        
    elseif isa(pow,'ADvar') == false
        adObj.der = pow .* adObj.val .^ (pow-1) .* adObj.der;
        adObj.val = adObj.val .^ pow;
        out = adObj;

    else
        out = exp(pow * log(adObj)); % calls ADvar's log()
    end
end