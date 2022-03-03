function J = spatialJac(varargin)
% SPATIALJAC(gst0, [{Y1,var1},{Y2,var2},...,{Yn,varn}]) computes the
%   Spatial Jacobian J of a serial manipulator.
%
%   Input
%       gst0:       initial configuration
%
%       varargin:   expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] the same
%                   number of 1x2 cells as the number of joints of the
%                   robot.
%                   Yn : n-th joint's unit twist
%                   varn : n-th joint variable
%

    [~, n] = size(varargin);
    
    g = expTw(varargin{1}{1}, varargin{1}{2});
    J = zeros(6, n);
    J(:, 1) = varargin{1}{1};

    for i = 2:n
        J(:, i) = adjoint(g) * varargin{i}{1};
        g = g * expTw(varargin{i}{1}, varargin{i}{2});
    end

end