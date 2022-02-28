function Js = SpJac(screwObj, varargin)
% SPJAC(screwObj, gst0, [{Y1,var1},{Y2,var2},...,{Yn,varn}]) computes the
%   Spatial Jacobian JS of a serial manipulator.
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
    
    Js = varargin{1}{1};

    if n == 1
        return
    end
    
    g = screwObj.expTw(varargin{1}{1}, varargin{1}{2});
    Js = [Js, zeros(6, n-1)];
    
    for i = 2:1:n
        Ytilde = screwObj.adjoint(g) * varargin{i}{1};
        Js = [Js(:, 1:i-1), Ytilde];
        g = g * screwObj.expTw(varargin{i}{1}, varargin{i}{2});
    end

end