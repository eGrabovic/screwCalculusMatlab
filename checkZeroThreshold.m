function iszero = checkZeroThreshold(tocheck, user_threshold)
% CHECKZEROTHRESHOLD(tocheck, user_threshold) checks if TOCHECK has all
%   elements that are "almost zero", within a tolerance, and returns a 
%   boolean ISZERO accordingly.
%   The tolerance can be specified with the argument USER_THRESHOLD or 
%   it's given a default value of 10^(-8).

    if exist('user_threshold', 'var')
        threshold = user_threshold;
    else
        threshold = 1e-8;
    end

    iszero = all(abs(tocheck) <= threshold, 'all');
end