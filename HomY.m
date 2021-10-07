function mat = HomY(theta, pos)
% TODO: if theta = 0 Rot=Id
    mat = [RotY(theta), [pos(1); pos(2); pos(3)]; ...
          [0, 0, 0],             1             ];
end