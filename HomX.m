function mat = HomX(theta, pos)
% HOMX(theta, pos) returns the homogeneous matrix generated through a
%   rotation of angle THETA around the X-axis and a translation POS

% TODO: if theta = 0 Rot=Id
    mat = [RotX(theta), [pos(1); pos(2); pos(3)]; ...
          [0, 0, 0],            1              ];
end