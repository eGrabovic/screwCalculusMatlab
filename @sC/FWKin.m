function gst = FWKin(screwObj, gst0,varargin)
% FWKIN(screwObj, gst0, varargin) computes the serial cinematic of a robot 
%   through Global P.O.E. parametrization.
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

    gst = sC.expTw(varargin{1}{1}, varargin{1}{2});

    if n > 1
        for i = 2:1:n
            gst = gst * sC.expTw(varargin{i}{1}, varargin{i}{2});
        end
    end

    gst = gst * gst0;
end