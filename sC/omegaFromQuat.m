function out = omegaFromQuat(q, qdot)
% TODO

    q0 = q(1);
    qvec = [q(2);q(3);q(4)];
    
    out = pinv(0.5.*[-qvec.'; -hat(qvec) + eye(3).*q0])*qdot;
end