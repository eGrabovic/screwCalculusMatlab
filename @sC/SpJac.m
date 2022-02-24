function Js = SpJac(varargin)
%
% Js = SpJac({Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola lo Jacobiano Spatial di un seriale
%
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto;
% varn : n esima variable di giunto;
%

[~,n] = size(varargin);

if n == 1
    Js = varargin{1}{1};
    return
end
    Js = varargin{1}{1};

g = sC.expTw(varargin{1}{1},varargin{1}{2});
Js = [Js zeros(6,n-1)];

for i = 2: 1 : n
    
    Ytilde = sC.adjoint(g)*varargin{i}{1};
    Js = [Js(:,1:i-1) Ytilde];
    g = g*sC.expTw(varargin{i}{1},varargin{i}{2});
end

end