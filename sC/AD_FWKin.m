function gst = AD_FWKin(gst0,varargin)
            %
            % Gst = FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
            % Funzione che calcola la cinematica seriale tramite
            % parametrizzazione GLOBAL P.O.E + la derivata tramite tecniche di
            % Automatic Differentiation.
            %
            % INPUTs:
            %
            % gst0 : offset tra spatial e tool quando i giunti sono nelle
            % condizioni iniziali.
            %
            % varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
            % giunti del seriale;
            % Yn : n esimo twist unitario del n esimo giunto;
            % varn : n esima variable di giunto immessa come variable di tipo
            % ADvar (proprietà con valore e derivata);
            %
            [~,n] = size(varargin);
            if n == 1
                gst = expTw(varargin{1}{1},varargin{1}{2});
                gst = gst*gst0;
                return
            end
            gst = expTw(varargin{1}{1},varargin{1}{2});
            for i = 2 : 1 : n
                gst = gst*expTw(varargin{i}{1},varargin{i}{2});
            end
            gst = gst*gst0;
        end