function mat = QuatToMat(quat)
% QUATTOMAT(quat) returns the rotation matrix MAT associated with the unit 
%   quaternion QUAT

    scalar = quat(1);
    vectorial = quat(2:end);
    
    Im = eye(3);
    vect_hat = Skew(vectorial);

    mat = Im + 2*vect_hat * (scalar * Im + vect_hat);
end