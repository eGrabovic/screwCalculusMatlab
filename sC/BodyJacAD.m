function Jb = BodyJacAD(gst0,varargin)
%
% Jb = BodyJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola lo Jacobiano Body di un seriale
%
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto;
% varn : n esima variable di giunto;
%
% gst0: configurazione di riferimento iniziale
%
[~,n] = size(varargin);
g = ADvar(gst0, zeros(4));

%             if n == 1
%                 Ytilde = adjoint(hom_mat_inv(g))*varargin{1}{1};
%                 Jb = Ytilde;
%                 return
%             end

Jb = ADvar(zeros(6, n), zeros(6, n));

for i = n : -1 : 1
    Ytilde = adjoint(hom_mat_inv(g))*varargin{i}{1};
    Jb.val(:,i) = Ytilde.val;
    Jb.der(:,i) = Ytilde.der;
    g = expTw(varargin{i}{1},varargin{i}{2})*g;
end

end