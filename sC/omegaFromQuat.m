function out = omegaFromQuat(q, qdot)
% OMEGAFROMQUAT(q, qdot) computes the angular velocity OUT given the
%   quaternion Q and its derivative QDOT.

    q0 = q(1);
    qvec = [q(2); q(3); q(4)];
    
    out = pinv(0.5 .* [-qvec.'; -hat(qvec) + eye(3) .* q0]) * qdot;
end