function uT = unitTwist(jointType,axis,q)
    %
    % Crea un Twist unitario a partire dal tipo di giunto ('P' per
    % prismatico o 'R' per rotoidale), l'asse del giunto ed
    % eventualmente la distanza dall origine q se di tipo 'R'
    %
    if strcmpi(jointType,'P') == true
        uT = [axis;0;0;0];
    end
    
    if strcmpi(jointType,'R') == true
        uT = [-cross(axis,q);axis];
    end
    
end