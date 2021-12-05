function quat = EulToQuat(EU_angles, type, deg_bool)
    quat = MatToQuat(EulToMat(EU_angles,type, deg_bool));
end