function quat = MatToQuat(mat)
% Warning: these are the reverse of conventional argument order:
% change order once all is consistent to (angle, nhat).
    diag_sum = trace(mat);

    if diag_sum > 0
        s = sqrt(diag_sum + 1.0);
        q0 = s/2;
        s = 1/(2*s);
        q1 = (mat(3,2) - mat(2,3))*s;
        q2 = (mat(1,3) - mat(3,1))*s;
        q3 = (mat(2,1) - mat(1,2))*s;
    else
        if (mat(1,1) >= mat(2,2)) && (mat(1,1) >= mat(3,3))
            % i=0,  j = 1, k = 2 *)
            s = sqrt(mat(1,1) - mat(2,2) - mat(3,3) + 1.0);
            q1 = s/2; s = 1/(2*s);
            q0 = (mat(3,2) - mat(2,3))*s;
            q2= (mat(2,1) + mat(1,2))*s;
            q3 = (mat(1,3) + mat(3,1))*s;
        else
            if (mat(1,1) < mat(2,2)) && (mat(1,1) >= mat(3,3))
                % i=1,  j = 2, k = 0 *)
                s = sqrt(mat(2,2) - mat(3,3) - mat(1,1) + 1.0);
                q2 = s/2; 
                s = 1/(2*s);
                q0 = (mat(1,3) - mat(3,1))*s;
                q3= (mat(3,2) + mat(2,3))*s;
                q1 = (mat(2,1) + mat(1,2))*s;
            else % (mat(1,1) < mat(2,2)) && (mat(1,1) < mat(3,3))
                % i=2,  j = 0, k = 1 *)
                s = sqrt(mat(3,3) - mat(1,1) - mat(2,2) + 1.0);
                q3 = s/2; 
                s = 1/(2*s);
                q0 = (mat(2,1) - mat(1,2))*s;
                q1= (mat(1,3) + mat(3,1))*s;
                q2 = (mat(3,2) + mat(2,3))*s;
            end
        end
    end
    quat = normalize([q0, q1, q2, q3]);
end