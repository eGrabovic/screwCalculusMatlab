function inv_quat = QuatInv(quat)
    inv_quat = [quat(1), - quat(2), -quat(3), - quat(4)];
end