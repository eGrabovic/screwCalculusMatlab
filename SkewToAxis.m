function vet = SkewToAxis(mat)
    % TODO: implement skewQ and use it here
    if any(size(mat) ~= [3,3]) % and not skewQ
        disp("Argument has wrong dimensions!");
        return
    end

    vet = [mat(3,2), mat(1,3), mat(2,1)];
end