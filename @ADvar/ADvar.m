classdef ADvar
    % ADvar : dual variable for automatic differentiation
    % 
    % This version supports correctly only 1 variable derivative; check "ADvar_multigradientsOPT.m" for multiple gradient derivative (not worth for this project cause of forward application)
    
    properties
        
        val % value
        der % derivative
        
    end
    
    methods
        
        %% COSTRUTTO/INIZIALIZZATORE
        function obj = ADvar(a,b)
            %
            % Costrutto/inizializzatore
            %
            if nargin == 0
                
                obj.val = [];
                obj.der = [];
                
            elseif nargin == 1      % caso in cui la variabile è una costante
                
                obj.val = a;
                obj.der = 0;
%                 [r,c] = size(a);      % this may be useful but it requires testing
%                 obj.der = zeros(r,c);
                
            else                     % caso generico in cui a è il valore e b è la derivata
                
                obj.val = a;
                obj.der = b;
                
            end
        end
        
        %% Programmo il comportamento delle funzioni built-in di matlab quando hanno in input una ADvar
         % + alcune funzioni di Screw_Calc
        
         function out = eq(u,v)
             
             if ~isa(u,'ADvar')
                 
                 out = u == v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val == v;
                 return
                 
             end
             
             out = u.val == v.val;
             
         end
         
         function out = ne(u,v)
             
             if ~isa(u,'ADvar')
                 
                 out = u ~= v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val ~= v;
                 return
                 
             end
             
             out = u.val ~= v.val;
             
         end
         
         function out = gt(u,v)
             
             if ~isa(u,'ADvar')
                 
                 out = u > v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val > v;
                 return
                 
             end
             
             out = u.val > v.val;
             
         end
         
         function out = ge(u,v)
             
             if ~isa(u,'ADvar')
                 
                 out = u >= v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val >= v;
                 return
                 
             end
             
             out = u.val >= v.val;
             
         end
         
         function out = lt(u,v)
             
             if ~isa(u,'ADvar')
                 
                 out = u < v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val < v;
                 return
                 
             end
             
             out = u.val < v.val;
             
         end
         
         function out = le(u,v)
             
             
             if ~isa(u,'ADvar')
                 
                 out = u <= v.val;
                 return
                 
             end
             
             if ~isa(v,'ADvar')
                 
                 out = u.val <= v;
                 return
                 
             end
             
             out = u.val <= v.val;
             
         end
         
        function u = uminus(u)
            
            u.val = -u.val;
            u.der = -u.der;
            
        end
         
        function u = cross(u,v)
            
           %out = ADvar(cross(u.val,v.val),cross(u.der,v.val) + cross(u.val,v.der));
           u.der = cross(u.der,v.val) + cross(u.val,v.der);
           u.val = cross(u.val,v.val);
            
        end
        
        function out = minus(u,v)
            
            if ~isa(u,'ADvar')
                % out = ADvar(u - v.val, v.der);
                v.der = -v.der;
                v.val = u - v.val;
                out = v;
            elseif ~isa(v,'ADvar')
                % out = ADvar(u.val - v, u.der);
                u.der = u.der;
                u.val = u.val - v;
                out = u;
            else
                % out = ADvar(u.val - v.val,u.der - v.der);
                u.der = u.der - v.der;
                u.val = u.val - v.val;
                out = u;
            end
        end
        
        function out = plus(u,v)
            %
            %
            %
            if isa(u,'ADvar') == false
                % out = ADvar(u + v.val, v.der);
                v.der = v.der;
                v.val = u + v.val;
                out = v;
            elseif isa(v,'ADvar') == false
                % out = ADvar(v + u.val, u.der);
                u.der = u.der;
                u.val = v + u.val;
                out = u;
            else
                % out = ADvar(u.val + v.val,u.der + v.der);
                u.der = u.der + v.der;
                u.val = u.val + v.val;
                out = u;
            end
        end
        
        function out = mtimes(u,v)
            % overload dell'operatore '*'
            %
            %
            
            if isa(u,'ADvar') == false
                
%                 out = ADvar(u*v.val,u*v.der);
                v.der = u*v.der;
                v.val = u*v.val;
                out = v;
                
            elseif isa(v,'ADvar') == false
                
%                 out = ADvar(u.val*v,u.der*v);
                u.der = u.der*v;
                u.val = u.val*v;
                out = u;
                
            else
%                 out = ADvar(u.val*v.val,u.der*v.val + u.val*v.der);
                u.der = u.der*v.val + u.val*v.der;
                u.val = u.val*v.val;
                out = u;
                
            end
        end
        
        function out = times(u,v)
            % overload dell'operatore '.*'
            %
            %
            
            if ~isa(u,'ADvar')
                
                % out = ADvar(u.*v.val,u.*v.der);
                v.der = u.*v.der;
                v.val = u.*v.val;
                out = v;
                
            elseif ~isa(v,'ADvar')
                
                %out = ADvar(u.val.*v,u.der.*v);
                u.der = u.der.*v;
                u.val = u.val.*v;
                out = u;
                
            else
                % out = ADvar(u.val.*v.val,u.der.*v.val + u.val.*v.der);
                u.der = u.der.*v.val + u.val.*v.der;
                u.val = u.val.*v.val;
                out = u;
                
            end
        end
        
        function out = mrdivide(u,v)
            
            if isa(u,'ADvar') == false
                
                % out = ADvar(u/v.val, -u.*v.der/v.val^2);
                v.der = -u.*v.der/v.val^2;
                v.val = u/v.val;
                out = v;
                
            elseif isa(v,'ADvar') == false
                
                % out = ADvar(u.val/v,u.der/v);
                u.der = u.der/v;
                u.val = u.val/v;
                out = u;
                
            else
                % out = ADvar(u.val/v.val,(u.der*v.val - u.val*v.der)/v.val^2);
                u.der = (u.der*v.val - u.val*v.der)/v.val^2;
                u.val = u.val/v.val;
                out = u;
                
            end
        end
        
        function out = rdivide(u,v)
            
            if isa(u,'ADvar') == false
                
                % out = ADvar(u./v.val, -u.*v.der./v.^2);
                v.der = -u.*v.der./v.^2;
                v.val = u./v.val;
                out = v;
                
            elseif isa(v,'ADvar') == false
                
                % out = ADvar(u.val./v,u.der./v);
                u.der = u.der./v;
                u.val = u.val./v;
                out = u;
                
            else
                % out = ADvar(u.val./v.val,(u.der.*v.val - u.val.*v.der)./v.val.^2);
                u.der = (u.der.*v.val - u.val.*v.der)./(v.val.^2);
                u.val = u.val./v.val;
                out = u;
                
            end
        end
        
        function out = mpower(u,v)
            %
            %
            %
            if isa(u,'ADvar') == false
                % out = ADvar(u^v.val,u^v.val*log(u)*v.der);
                v.der = u^v.val*log(u)*v.der;
                v.val = u^v.val;
                out = v;
                
            elseif isa(v,'ADvar') == false
                % out = ADvar(u.val^v, v*u.val^(v-1)*u.der);
                u.der = v*u.val^(v-1)*u.der;
                u.val = u.val^v;
                out = u;
                
            else
                out = exp(v*log(u));
            end
        end
        
        function out = power(u,v)
            %
            %
            %
            if isa(u,'ADvar') == false
                % out = ADvar(u.^v.val,u.^v.val*log(u)*v.der);
                v.der = u.^v.val.*log(u).*v.der;
                v.val = u.^v.val;
                out = v;
                
            elseif isa(v,'ADvar') == false
                % out = ADvar(u.val.^v, v*u.val.^(v-1)*u.der);
                u.der = v.*u.val.^(v-1).*u.der;
                u.val = u.val.^v;
                out = u;
            else
                out = exp(v*log(u));
            end
        end
        
        function out = double(u)
            %
            %
            %
            out = [u.val, u.der];
            
        end
        
        function u = sqrt(u)
            %
            % overload funzione radcie
            %
            
            %u = ADvar(sqrt(u.val),1./(2*sqrt(u.val))*u.der);
            u.der = 1./(2*sqrt(u.val)).*u.der;
            u.val = sqrt(u.val);
            
        end
        
        function u = sin(u)
            % overload funzone seno
            %
            %
            
            %out = ADvar(sin(u.val),cos(u.val).*u.der);
            u.der = cos(u.val).*u.der;
            u.val = sin(u.val);
            
        end
        
        function u = asin(u)
            %
            %
            %
            
            %out = ADvar(asin(u.val),1/sqrt(1-u.val*u.val)*u.der);
            u.der = 1./sqrt(1-u.val.*u.val).*u.der;
            u.val = asin(u.val);
            
        end
        
        function u = cos(u)
            % overload funzione coseno
            %
            %
            
            % out = ADvar(cos(u.val),-sin(u.val).*u.der);
            u.der = -sin(u.val).*u.der;
            u.val = cos(u.val);
            
        end
        
        function u = acos(u)
            %
            %
            %
            
            %out = ADvar(acos(u.val),-1./sqrt(1-u.val*u.val)*u.der);
            u.der = -1./sqrt(1-u.val*u.val)*u.der;
            u.val = acos(u.val);
            
        end
        
        function u = tan(u)
            % tan overload
            %
            %
            
            % out = sin(u)./cos(u);
            u.der = 2./(cos(2.*u.val) + 1);
            u.val = tan(u.val);
            
        end
        
        function Y = atan2(Y,X)
            %
            %
            %
            
            %out = ADvar(atan2(Y.val,X.val), X.val*Y.der/(X.val.^2 + Y.val.^2) - Y.val*X.der/(X.val.^2 + Y.val.^2));
            Y.der = X.val*Y.der/(X.val.^2 + Y.val.^2) - Y.val*X.der/(X.val.^2 + Y.val.^2);
            Y.val = atan2(Y.val,X.val);
            
            
        end
        
        function u = atan(u)
            %
            %
            %
            
            %out = ADvar(atan2(Y.val,X.val), X.val*Y.der/(X.val.^2 + Y.val.^2) - Y.val*X.der/(X.val.^2 + Y.val.^2));
            u.der = u.der./(u.val.^2 + 1);
            u.val = atan(u.val);
            
            
        end
        
        function u = log(u)
            % overload log
            %
            %
            
            % out = ADvar(log(u.val),u.der/u.val);
            u.der = u.der./u.val;
            u.val = log(u.val);
            
        end
        
        function u = exp(u)
            % exp overload
            %
            %
            
%             out = ADvar(exp(u.val),u.der.*exp(u.val));
            u.der = u.der.*exp(u.val);
            u.val = exp(u.val);
            
        end
        
        function u = abs(u)
            % abs overload
            %
            %
            
            %out = ADvar(abs(u.val),u.der*sign(u));
            u.der = u.der*sign(u.val);
            u.val = abs(u.val);
            
        end
        
        function out = getVal(u)
            %
            %
            %
            out = u.val;
            
        end
        
        function out = getDer(u)
            %
            %
            %
            out = u.der;
            
        end
        
        function out = horzcat(varargin)
            
            if ~isa(varargin{2},'ADvar')
                varargin{2} = ADvar(varargin{2}, zeros(size(varargin{2})));
            elseif ~isa(varargin{1},'ADvar')
                varargin{1} = ADvar(varargin{1}, zeros(size(varargin{1})));
            end
            obj = cat(2, varargin{:});
            as  = horzcat(obj.val);
            bs  = horzcat(obj.der);
            out = ADvar( as, bs);
            
        end
         
        function out = vertcat(varargin)

            for k = 1:length(varargin) 
                if ~isa(varargin{k}, 'ADvar')
                    varargin{k} = ADvar(varargin{k}, zeros(size(varargin{k})));
                end
            end
            obj = cat(1, varargin{:});
            as  = vertcat(obj.val);
            bs  = vertcat(obj.der);

            out = ADvar( as, bs );
            
        end

        
%         function out = cat(dim, varargin)
%             
%             obj = [varargin{:}];
%             
%             a = cat(dim, obj.val);
%             b = cat(dim, obj.der);
%             out = ADvar(a,b);
%             
%         end
        
        function A = transpose(A)
            %
            %
            %
            %out = ADvar(A.val.',A.der.');
            A.der = A.der.';
            A.val = A.val.';
            
        end
        
        function u = subsref(u,s)
            switch s(1).type
                case '.'
                    u = builtin('subsref',u,s);
                case '()'
                    u.der = subsref(u.der, s);
                    u.val = subsref(u.val, s);
            end
        end
        
        function u = subsasgn(u,s,val)
            switch s(1).type
                case '.'
                    u = builtin('subsasgn',u,s,val);
                case '()'
                    if isa(val, 'ADvar')
                        u.der = subsasgn(u.der, s, val.der);
                        u.val = subsasgn(u.val, s, val.val);
                    else
                        u.der = subsasgn(u.der, s, 0);
                        u.val = subsasgn(u.val, s, val);
                    end
            end
        end
        
        function u = sum(u, dim)
            
            u.val = sum(u.val, dim);
            u.der = sum(u.der, dim);

        end
        
    end
end