function Y = Regressor(DHtable, q, qp, v, vp, t, g0, k)
% REGRESSOR(DHtable, q, qp, v, vp, t, g0, k) directly computes the 
%   regressor for the supplied: 
%   - DHTABLE (Denavit-Hartemberg table)
%   - Q (configuration)
%   - QP (velocity)
%   - V (reference velocity)
%   - VP (reference acceleration)
%   - T (independent var)
%   - G0 (components of gravity in {0})
%   - K (contribute for the k-th link). 
%   If k is omitted, the complete regressor is returned

    R0k  = RigidOrientation(DHFKine(DHtable, k));
    R0kT = R0k';

    Jk   =  DHJacob0Dyn(DHtable, k);
    Jvk  =  Jk(1:3, :);
    Jok  =  Jk(4:6, :);
    
    JvkT = Jvk';
    JokT = Jok';
    JvTJv = JvkT * Jvk;
    
    Jvkp  = TensorDerivative(Jvk, q) * qp;
    JvkpT = Jvkp';
    
    Jokp  = TensorDerivative(Jok, q) * qp;
    JokpT = Jokp'; 
		    
    % termini da d   dT1
    %            -   -  
    %            dt  dqp         con T1 = v*B(q)*qp
    
    X0p = ((TensorDerivative(JvTJv, q) * qp) * v  + JvTJv * vp)';
    
    
    % X1p = (JvkpT * Hat(Jok * v) + JvkT * (Hat(Jokp * v) + Hat(Jok * vp)) - ...
    % JokpT * Hat(Jvk * v)  - JokT * (Hat(Jvkp * v) + Hat(Jvk * vp)) ) * R0k +...
    % (JvkT * Hat(Jok * v) - JokT * Hat(Jvk * v)) * (diff(R0k, t));
    
    T1 = sparse([2 3 3],[3 2 3],[-1 1 0]);
    T2 = sparse([1 3 3],[3 1 3],[1 1 0]);
    T3 = sparse([1 2 3],[2 1 3],[-1 1 0]);
    
    TT = [T1; T2; T3];
    
    X1a  = JokT * R0k * TT * R0kT * Jvk;
    
    X1ap = [];
    for i = [1,3]
        X1ap = [X1ap; TensorDerivative(X1a(:,i,:), q) * qp * v];
    end
    
    X1b  = -JvkT * R0k * TT * R0kT * Jok;  
    
    X1bp = [];
    for i = [1,3]
        X1bp = [X1bp; TensorDerivative(X1b(:,i,:), q) * qp * v];
    end
    
    X1p = (X1a + X1b) * vp + X1ap + X1bp;
    
    E1 = sparse([1 3],[1 3],[1 0]);
    E2 = sparse([1 2 3],[2 1 3],[-1 -1 0]);
    E3 = sparse([1 3 3],[3 1 3],[-1 -1 0]);
    E4 = sparse([2 3],[2 3],[1 0]);
    E5 = sparse([2 3 3],[3 2 3],[-1 -1 0]);
    E6 = sparse(3, 3, 1);
    EE = [E1; E2; E3; E4; E5; E6];
    
    par = JokT * R0k * EE * R0kT * Jok;
    
    X2p = [];
    for i = [1,6]
        X2p = [X2p; TensorDerivative(par(:,i,:), q) * qp * v + par * vp];
    end

    % termini da dT2
    %            -
    %            dq           con T2 = (1/2)*v*B(q)*qp
    
    W0 = (1/2) * diff( v * JvkT * Jvk * qp, q);
    
    W1 = (1/2) * (diff((v*JvkT * Hat(Jok*qp) - v*JokT * Hat(Jvk*qp)) * R0k, q))';
    
    W2 = (1/2) * (diff(v * JokT * R0k * EE * R0k' * Jok * qp, q))';
    
    % termini da d U
    %            -
    %            d q
    
    Z0 = -JvkT * g0;
    
    Z1 = -(diff(g0*R0k, q))'; 
    
    Y0 = X0p - W0  + Z0;
    Y1 = X1p - W1  + Z1;
    Y2 = X2p - W2;
    
    Y  = [Y0, Y1, Y2];
end