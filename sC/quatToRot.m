function out = quatToRot(Quat)
% TODO

q0 = Quat(1);
q = zeros(3,1, class(Quat(1)));
q(1) = Quat(2);
q(2) = Quat(3);
q(3) = Quat(4);

out = eye(3) + 2*hat(q)*q0 + 2*hat(q)*(q0.*eye(3) + hat(q));


end