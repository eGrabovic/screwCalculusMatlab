function E = expTw(unitTwist,th)
% EXPTW(unitTwist, th) computes the exponential E of a twist UNITTWIST w.r.t.
%   the joint variable TH.
    
    if isa(th, 'ADvar')
       E = expTw(unitTwist, th.val);
       E = ADvar(E, hat(unitTwist)*E.*th.der);
       return
    end
    
    if all(isnumeric(unitTwist)) == true
        if abs(unitTwist(4)) <= 1e-6 && abs(unitTwist(5)) <= 1e-6 && abs(unitTwist(6)) <= 1e-6 % prismatic joint
            
            axis = [unitTwist(1);unitTwist(2);unitTwist(3)];
            R = eye(3);
            d = axis*th;
            
        else % revolute joint
            
            axis = [unitTwist(4);unitTwist(5);unitTwist(6)];
            v = [unitTwist(1);unitTwist(2);unitTwist(3)];
            q = -cross(v,axis);
            R = expSkew(axis,th);
            d =(eye(3) - R)*q;
            
        end
%                 E = hom_mat(R,d);
        E = [R d; 0 0 0 1];

    end
end