function X = adjoint(screwObj, rotm,disp)
% ADJOINT(Gst) computes the adjoint transformation X of a 4x4 homogeneous 
%   matrix ROTM or of a 3x3 rotation matrix ROTM and a displacement vector 
%   DISP.
%   X is a 6x6 matrix that can be used as a congruent twists.

    % Only a homogeneous matrix
    if nargin == 1 % TODO: does this increase because of the class object?
        if all(size(rotm) ~= [4,4])
            error('Gst homogeneous matrix must be a 4x4 matrix')
        end

        R = rotm(1:3,1:3);
        d = rotm(1:3,4);

        X = [[R,            screwObj.hat(d) * R];
             [zeros(3, 3),  R]];

    else % Rotation matrix + displacement
        if all(size(rotm) ~= [3,3]) || all(size(disp) ~= [3,1])
            error(['Rotation matrix must be a 3x3 matrix and the ' ...
                'displacement vector a 3x1 array'])
        end

        R = rotm;
        d = disp;

        X = [[R,            screwObj.hat(d) * R];
             [zeros(3, 3),  R]];
    end
end