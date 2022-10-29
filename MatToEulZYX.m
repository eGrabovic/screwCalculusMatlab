function eul_angles = MatToEulZYX(Rotmat, deg_bool)
% MATTOEULZYZ(Rotmat) gathers the Euler angles EUL_ANGLES from the supplied Rotation
%   matrix ROTMAT.
%   The angle theta is considered to be in [-pi/2; pi/2].
%   If DEG_BOOL is specified as true, output angles will be given in
%   degrees.
%
%   Examples
%       MatToEulZYX([[1    0         0]
%                    [0    0.7071   -0.7071]
%                    [0    0.7071    0.7071]])
%   --> 0   0   0.7854
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
    
	% Check the singularity of representation ZYX
	thetaisplushalfPi = abs(Rotmat(3,1) + 1) < 10^(-10);
	thetaisminushalfPi = abs(Rotmat(3,1) - 1) < 10^(-10);     
	        
	% In cases of singularity we set arbitrarily --> psi = 0
    if thetaisplushalfPi
	     phi = atan2(Rotmat(2,3), Rotmat(1,3));
	     theta = pi/2;
	     psi = 0;
    end
	  
    if thetaisminushalfPi
	     phi = atan2(-Rotmat(2,3), -Rotmat(1,3));
	     theta = -pi/2;
	     psi = 0;
    end
	   
    if ~(thetaisplushalfPi || thetaisminushalfPi)
	     phi = atan2(Rotmat(2,1), Rotmat(1,1));
	     theta = atan2(-Rotmat(3,1), sqrt(Rotmat(3,2)^2 + Rotmat(3,3)^2));
	     psi = atan2(Rotmat(3,2), Rotmat(3,3));
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
