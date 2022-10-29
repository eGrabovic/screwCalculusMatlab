function H = JointCenterRange(q, q_min, q_max)

    assert(isvector(q), "q is not a vector");
    assert(isvector(q_min), "q_min is not a vector");
    assert(isvector(q_max), "q_max is not a vector");

    assert(all(size(q) == size(q_min)) || all(size(q) == circshift(size(q_min),[0 1])), ...
      "Minimum limit array and supplied joint values have inconsistent " + ...
      "size: q: %s q_min: %s", mat2str(size(q)), mat2str(size(q_min)));
    assert(all(size(q) == size(q_max)) || all(size(q) == circshift(size(q_max),[0 1])), ...
      "Mamimum limit array and supplied joint values have inconsistent " + ...
      "size: q: %s q_max: %s", mat2str(size(q)), mat2str(size(q_max)));

    % Since assert did not trigger, arrays have the same size or are
    % transposed
    if all(size(q) == circshift(size(q_min),[0 1])) % q: [x y] & q_min: [y x]
        qm = q_min';
    else
        qm = q_min;
    end

    if all(size(q) == circshift(size(q_max),[0 1])) % q: [x y] & q_max: [y x]
        qM = q_max';
    else
        qM = q_max;
    end

    n = max(size(q));
    [a, b] = size(q);

    q_half = (qm + qM) / 2;
    q_diff = qM - qm;

    H = reshape(1/(2*n) * sum((q - q_half).^2 / q_diff.^2), [a b]);
end