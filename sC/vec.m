function X = vec(Tw)
%
% transforms a hat form twist into the vec form twist: R_4x4 -> R_6
%
%

if isa(Tw, 'ADvar')
    X = ADvar(vec(Tw.val), vec(Tw.der)); % recursively apply the function to the value and derivative parts separately
    return
end

if size(Tw, 1) == 3
    
    X(1) = Tw(3,2);
    X(2) = Tw(1,3);
    X(3) = Tw(2,1);
    X = X.';
    return
end

w_hat = Tw(1:3,1:3);
v = Tw(1:3,4);
X = [v; vec(w_hat)];

end