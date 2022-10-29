function quat = MatToQuat(mat)
% MATTOQUAT(R) returns the unit quaternion QUAT associated with the 
%   rotation matrix MAT

% WARNING: these are the reverse of conventional argument order;
% change order once all is consistent to (angle, nhat). TODO?
% MATLAB uses [w; x,y,z]
% here: 
    diag_sum = trace(mat);

    if diag_sum > 0
        s = sqrt(diag_sum + 1.0);
        qw = s/2;
        s = 1/(2*s);
        qx = (mat(3,2) - mat(2,3))*s;
        qy = (mat(1,3) - mat(3,1))*s;
        qz = (mat(2,1) - mat(1,2))*s;

    else
        if (mat(1,1) >= mat(2,2)) && (mat(1,1) >= mat(3,3))
            s = sqrt(mat(1,1) - mat(2,2) - mat(3,3) + 1.0);
            qx = s/2; 
            s = 1/(2*s);
            qw = (mat(3,2) - mat(2,3))*s;
            qy= (mat(2,1) + mat(1,2))*s;
            qz = (mat(1,3) + mat(3,1))*s;

        else
            if (mat(1,1) < mat(2,2)) && (mat(1,1) >= mat(3,3))
                s = sqrt(mat(2,2) - mat(3,3) - mat(1,1) + 1.0);
                qy = s/2; 
                s = 1/(2*s);
                qw = (mat(1,3) - mat(3,1))*s;
                qz= (mat(3,2) + mat(2,3))*s;
                qx = (mat(2,1) + mat(1,2))*s;

            else % (mat(1,1) < mat(2,2)) && (mat(1,1) < mat(3,3))
                s = sqrt(mat(3,3) - mat(1,1) - mat(2,2) + 1.0);
                qz = s/2; 
                s = 1/(2*s);
                qw = (mat(2,1) - mat(1,2))*s;
                qx= (mat(1,3) + mat(3,1))*s;
                qy = (mat(3,2) + mat(2,3))*s;
            end
        end
    end
    quat = normalize([qw, qx, qy, qz]);
end