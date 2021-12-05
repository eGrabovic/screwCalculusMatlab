function or_err = MatToVectParam(Rotmat)
% MATTOVECT(R_ACT) returns the orientation error OR_ERR = r*sin(theta), where
%   r and theta are the components of the axis-angle parametrization.
%   OR_EE represents also the axis form of the skew-symmetrix matrix 
%   R_SS = (1/2)*(R - R^T).

    % TODO: how to check if it's a rotation matrix?
    assert(ismatrix(Rotmat) && all(size(Rotmat) == [3,3]),...
        "Wrong matrix dimensions");

	[axis, theta] = RotationParam(Rotmat);
	or_err = axis * sin(theta);
end 