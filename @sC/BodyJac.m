function Jb = BodyJac(screwObj, gst0, varargin)
% BODYJAC(screwObj, gst0, varargin) computes the body jacobian JB of a 
%   serial robot.
%  
%   Input
%       gst0:       initial configuration homogeneous matrix
%
%       varargin :  expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] with the same
%                   number of 1x2 cells as the number of joints of the robot.
%                   - Yn: n-th joint's unit twist
%                   - varn: n-th joint variable

    [~, n] = size(varargin);
    g = gst0;
    
    if n == 1
        Jb = screwObj.adjoint(screwObj.hom_mat_inv(g)) * varargin{1}{1};
        return
    end
    
    Jb = zeros(6, n);
    
    for i = n:-1:1
        Jb(:, i) = screwObj.adjoint(screwObj.hom_mat_inv(g)) * varargin{i}{1};
        g = screwObj.expTw(varargin{i}{1}, varargin{i}{2}) * g;
    end

end