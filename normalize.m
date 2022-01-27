function qnorm = normalize(q)
% NORMALIZE(q) ...
% TODO: check and fill

    norm = q*q';
    eps = 10^(-14); % tolerance

    if ~isnan(norm)
        if norm<eps
            qnorm = q;
        else
            qnorm = q/sqrt(norm);
        end
    else
        qnorm = q/sqrt(norm); % TODO: why?
    end
end