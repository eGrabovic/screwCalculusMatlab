%% script for spatial jacobian derivative computation
clc
clear all
% numJoints
n = 6;

% remember linear contibution of revolute joint = -cross(w, q)
% t = type: 0 -> prismatic; 1 -> revolute
w(:, 1) = [0;0;1]; q(:, 1) = [0;0;  0]; t(1) = 1;
w(:, 2) = [0;0;1]; q(:, 2) = [0;0;  0]; t(2) = 0;
w(:, 3) = [0;1;0]; q(:, 3) = [0;0;  2]; t(3) = 0;
w(:, 4) = [0;1;0]; q(:, 4) = [0;0;  2]; t(4) = 1;
w(:, 5) = [1;0;0]; q(:, 5) = [0;4;  2]; t(5) = 1;
w(:, 6) = [0;0;1]; q(:, 6) = [0;4.5;0]; t(6) = 1;

Y = zeros(6, n);
for i = 1:n % unitary twists computation
    if t(i) == 1 % revolute
        Y(:, i) = [-cross(w(:, i), q(:, i)); w(:, i)];
    elseif t(i) == 0 % prismatic
        Y(:, i) = [w(:, i); [0; 0; 0]];
    end
end

%% spatial Jacobian computation
% gst0 = hom_mat(eye(3), [0; 4; 2]);
% Jac_spatial = @(th) spatialJac(gst0,...
%     {Y(:, 1), th(1)},...
%     {Y(:, 2), th(2)},...
%     {Y(:, 3), th(3)},...
%     {Y(:, 4), th(4)},...
%     {Y(:, 5), th(5)},...
%     {Y(:, 6), th(6)});
% 
% q = [1,1,1,1,1,1];
% qdot = ones(1,6);
% disp(" jacobian computation with algebraic algorithmic evaluation ")
% JacDot = @(th, thdot) spatialJacDerivative(Jac_spatial(th), thdot);
% tic
% res = JacDot(q, qdot);
% toc
% res
% disp(" jacobian computation with finite difference ")
% 
% tic
% delta = 0.00001;
% deltaVec = zeros(1, 6);
% Jdot = zeros(6, n);
% for i = 1:n
%     dir = deltaVec;
%     dir(i) = delta;
%     Jdot = Jdot + (Jac_spatial(q + dir) - Jac_spatial(q))./delta.*qdot(i);
%     
% end
% toc
% Jdot

%% body Jacobian computation
gst0 = hom_mat(eye(3), [0; 4; 2]);
Jac_body = @(th) BodyJac(gst0,...
    {Y(:, 1), th(1)},...
    {Y(:, 2), th(2)},...
    {Y(:, 3), th(3)},...
    {Y(:, 4), th(4)},...
    {Y(:, 5), th(5)},...
    {Y(:, 6), th(6)});

q = [1,1,1,1,1,1];
qdot = ones(1,6);
disp(" body jacobian computation with algebraic algorithmic evaluation ")
JacDot = @(th, thdot) bodyJacDerivative(Jac_body(th), thdot);
tic
for s = 1:10000
    res = JacDot(q, qdot);
end
toc
res
disp(" body jacobian computation with finite difference ")

tic
delta = 0.00001;
deltaVec = zeros(1, 6);

for s = 1:10000
    Jdot = zeros(6, n);
    for i = 1:n
        dir = deltaVec;
        dir(i) = delta;
        Jdot = Jdot + (Jac_body(q + dir) - Jac_body(q))./delta.*qdot(i);
        
    end
end
toc
Jdot

Jb = @(th1, th2, th3, th4, th5, th6) BodyJacAD(gst0,...
    {Y(:, 1), th1},...
    {Y(:, 2), th2},...
    {Y(:, 3), th3},...
    {Y(:, 4), th4},...
    {Y(:, 5), th5},...
    {Y(:, 6), th6});

tic
for s = 1:10000
    res = Jb(ADvar(q(1),1),...
             ADvar(q(2),1),...
             ADvar(q(3),1),...
             ADvar(q(4),1),...
             ADvar(q(5),1),...
             ADvar(q(6),1));
end
res
toc
