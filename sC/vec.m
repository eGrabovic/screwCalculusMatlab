function X = vec(Tw)
% VEC(Tw) transforms a hat-form twist TW into the vec-form twist X.
%   Allows for ADvar class TW argument.

    % Manage instances of ADvar
    if isa(Tw, 'ADvar')
        % Recursively apply the function to the value and derivative parts 
        % separately
        X = ADvar(vec(Tw.val), vec(Tw.der)); 
        return
    end
    
    % Case of a 3x3 matrix
    if all(size(Tw) == [3, 3])
        X(1) = Tw(3,2);
        X(2) = Tw(1,3);
        X(3) = Tw(2,1);
        X = X.'; % To return a column vector
        return
    end
    
    % Final case: 4x4 homogeneous matrix

    assert(all(size(mat) == [4, 4]), ...
        "Supplied argument must be an ADvar instance or a 3x3 or 4x4 matrix");

    w_hat = Tw(1:3,1:3);
    v = Tw(1:3,4);

    % Return the whole twist
    X = [v; vec(w_hat)]; 

end