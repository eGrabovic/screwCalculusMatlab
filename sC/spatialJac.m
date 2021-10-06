function J = spatialJac(gst0, varargin)
%
% Jb = spatialJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola lo Jacobiano spatial di un seriale
%
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto;
% varn : n esima variable di configurazione del giunto;
%
% gst0: configurazione di riferimento iniziale
%
[~, n] = size(varargin);

g = expTw(varargin{1}{1}, varargin{1}{2});
J = zeros(6, n);
J(:, 1) = varargin{1}{1};
for i = 2:n
    J(:, i) = adjoint(g)*varargin{i}{1};
    g = g*expTw(varargin{i}{1}, varargin{i}{2});
end

end