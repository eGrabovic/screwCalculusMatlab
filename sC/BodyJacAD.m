function Jb = BodyJacAD(gst0,varargin)
% BODYJAC(gst0,varargin) computes the body jacobian JB of a serial robot.
%  
%   Input
%       gst0:       initial configuration homogeneous matrix
%
%       varargin :  expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] the same
%                   number of 1x2 cells as the number of joints of the robot.
%                   Yn: n-th joint's unit twist
%                   varn: n-th joint variable as a ADvar variable
%

[~,n] = size(varargin);
g = ADvar(gst0, zeros(4));

%             if n == 1
%                 Ytilde = adjoint(hom_mat_inv(g))*varargin{1}{1};
%                 Jb = Ytilde;
%                 return
%             end

Jb = ADvar(zeros(6, n), zeros(6, n));

for i = n : -1 : 1
    Ytilde = adjoint(hom_mat_inv(g))*varargin{i}{1};
    Jb.val(:,i) = Ytilde.val;
    Jb.der(:,i) = Ytilde.der;
    g = expTw(varargin{i}{1},varargin{i}{2})*g;
end

end