function gst = AD_FWKin(screwObj, gst0,varargin)
% AD_FWKIN(screwObj, gst0, varargin) computes the serial cinematic of a robot through 
%   Global P.O.E. parametrization and its derivate through Automatic 
%   Differentiation technics.
%
%   Input
%       gst0:       spatial-tool offset at initial conditions
%
%       varargin:   expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] the same
%                   number of 1x2 cells as the number of joints of the
%                   robot.
%                   Yn : n-th joint's unit twist
%                   varn : n-th joint variable as an ADvar instance
%
%   Output
%

    [~, n] = size(varargin);

    gst = ADexpTw(varargin{1}{1}, varargin{1}{2});
    
    if n > 1
        for i = 2 : 1 : n
            gst = gst * ADexpTw(varargin{i}{1}, varargin{i}{2});
        end
    end
    
    gst = gst * gst0;
end