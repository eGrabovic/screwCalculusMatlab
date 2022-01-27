function esp = TwistExp(xi, theta)
% TWISTEXPE(xi, Theta) gives the matrix exponential ESP of a twist XI.
%   Default value of THETA is 1.

% Take the exponential of a twist
% ! This only works in dimension 3 for now !

    w = xitow(xi);
    v = xitov(xi);
    
    % Make sure that we got a real twist
    assert(all(size(w) ~= 0) && all(size(v) ~= 0), ...
           "The provided vector is not a real twist");
    
    % Use the exponential formula from MLS
    if  all(w == [0,0,0]) || all(w == [0; 0; 0])
        Rot = eye(3);
        pos = v * theta;
    else
        ws = Skew(w);
        Rot = SkewExp(ws, theta);
        pos = (eye(3) - Rot) * (ws * v) + w * (w * v) * theta;
    end

    esp = RPToHomogeneous(Rot, pos);
end