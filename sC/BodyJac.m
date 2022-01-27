function Jb = BodyJac(gst0,varargin)
% BODYJAC(gst0,varargin) computes the body jacobian JB of a serial robot.
%  
%   Input
%       gst0:       initial configuration homogeneous matrix
%
%       varargin :  expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] the same
%                   number of 1x2 cells as the number of joints of the robot.
%                   Yn: n-th joint's unit twist
%                   varn: n-th joint variable
%

[~, n] = size(varargin);
cl = class(varargin{2}{2});
g = gst0;
Jb = zeros(6, n, cl);

for i = n : -1 : 1
    
    g = expTw(varargin{i}{1},varargin{i}{2})*g;  %g_i+1,n
    Jb(:, i) = adjoint(hom_mat_inv(g))*varargin{i}{1}; % adjoint(g_i,n^-1)*Y_i
    
end

end