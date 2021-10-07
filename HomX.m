function mat = HomX(theta, pos)
    mat = [RotX(theta), [pos(1); pos(2); pos(3)]; ...
          [0, 0, 0],            1              ];
end