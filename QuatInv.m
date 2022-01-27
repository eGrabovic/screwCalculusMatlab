function inv_quat = QuatInv(quat)
% QUATINV(quat) returns the conjugate (which is also the inverse) INV_QUAT
%   of a unit quaternion QUAT

    inv_quat = [quat(1), - quat(2), -quat(3), - quat(4)];
end