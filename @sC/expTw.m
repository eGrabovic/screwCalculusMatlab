function E = expTw(unitTwist,th)
    %
    % Calcolo dell esponenziale di un twist
    %
    % unitTwist = input twist unitario relativo al giunto
    %
    % th = variabile di giunto
    %
    if all(isnumeric(unitTwist)) == true
        if abs(unitTwist(4)) <= 1e-6 && abs(unitTwist(5)) <= 1e-6 && abs(unitTwist(6)) <= 1e-6
            
            axis = [unitTwist(1);unitTwist(2);unitTwist(3)];
            R = eye(3);
            d = axis*th;
            
        else
            
            axis = [unitTwist(4);unitTwist(5);unitTwist(6)];
            v = [unitTwist(1);unitTwist(2);unitTwist(3)];
            q = -cross(v,axis);
            R = sC.expSkew(axis,th);
            d =(eye(3) - R)*q;
            
        end
        E = sC.hom_mat(R,d);

    end
end