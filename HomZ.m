function mat = HomZ(theta, pos)
% HOMZ(theta, pos) returns the homogeneous matrix MAT generated through a
%   rotation of angle THETA around the Z-axis and a translation POS

% TODO: if theta = 0 Rot=Id
    mat = [RotZ(theta), [pos(1); pos(2); pos(3)]; ...
          [0, 0, 0],             1             ];
end