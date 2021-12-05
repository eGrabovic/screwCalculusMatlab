function eul_angles = MatToEulZYZ(Rotmat, deg_bool)
% MATTOEULZYZ(Rotmat) gathers the Euler angles EUL_ANGLES from the supplied Rotation
%   matrix ROTMAT.
%   The angle theta is considered to be in [0; pi].
%   If DEG_BOOL is specified as true, output angles will be given in
%   degrees.
%
%   Examples
%       MatToEulZYZ([[1    0         0]
%                    [0    0.7071   -0.7071]
%                    [0    0.7071    0.7071]])
%   --> -1.5708    0.7854    1.5708
%
%   Input
%       Rotmat:     3x3 matrix
%       deg_bool:   boolean
%   Output
%       eul_angles: 1x3 numerical array
%

    % TODO: how to check if it's a rotation matrix?
    assert(all(size(Rotmat) == [3,3]),...
           "Supplied matrix must be a 3x3 Rotation matrix");
	        
	% Check the singularity of representation ZYZ        
	thetaiszero = abs(Rotmat(3,3) - 1) < 10^(-10);
	thetaisPi = abs(Rotmat(3,3) + 1) < 10^(-10);     
	        
	% In cases of singularity we set arbitrarily --> psi = 0     
    if thetaiszero
	     phi = atan2(Rotmat(2,1) , Rotmat(1,1));
	     theta = 0;
	     psi = 0;
    end
      
    if thetaisPi
	     phi = atan2(Rotmat(2,1), Rotmat(1, 1));
	     theta = Pi;
	     psi = 0;
    end

    if ~(thetaiszero || thetaisPi)
	     phi = atan2(Rotmat(2,3), Rotmat(1,3));
	     theta = atan2(sqrt(Rotmat(1,3)^2 + Rotmat(2,3)^2 ), Rotmat(3,3));
	     psi = atan2(Rotmat(3,2), -Rotmat(3,1));
    end

    if (nargin == 2)
        if deg_bool
            eul_angles = [rad2deg(phi), rad2deg(theta), rad2deg(psi)];
        else
            eul_angles = [phi, theta, psi];
        end
    else
        eul_angles = [phi, theta, psi];
    end
	     
end