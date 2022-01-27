function gst = FWKin(gst0,varargin)
% AD_FWKIN(gst0,varargin) computes the serial cinematic through Global 
%   P.O.E. parametrization and its derivate through Automatic 
%   Differentiation technics.
%
%   Input
%       gst0:       spatial-tool offset at initial conditions
%
%       varargin:   expects [{Y1,var1},{Y2,var2},...,{Yn,varn}] the same
%                   number of 1x2 cells as the number of joints of the
%                   robot.
%                   Yn : n-th joint's unit twist unitario;
%                   varn : n-th joint variable as a ADvar instance
%
%   Output
%

n = length(varargin);
gst = expTw(varargin{1}{1},varargin{1}{2});
for i = 2 : n
    gst = gst*expTw(varargin{i}{1},varargin{i}{2});
end
gst = gst*gst0;
end