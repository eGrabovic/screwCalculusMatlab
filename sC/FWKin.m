function gst = FWKin(gst0,varargin)
%
% Gst = FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola la cinematica seriale tramite
% parametrizzazione GLOBAL P.O.E.
%
% INPUTs:
%
% gst0 : offset tra spatial e tool quando i giunti sono nelle
% condizioni iniziali.
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto;
% varn : n esima variable di giunto;
%
n = length(varargin);
gst = expTw(varargin{1}{1},varargin{1}{2});
for i = 2 : n
    gst = gst*expTw(varargin{i}{1},varargin{i}{2});
end
gst = gst*gst0;
end