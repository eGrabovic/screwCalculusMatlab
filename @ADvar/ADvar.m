classdef ADvar
% Dual variable for automatic differentiation
% 
% This version correctly supports only 1 variable derivative. 
% Check "ADvar_multigradientsOPT.m" for multiple gradient derivative 
% (not worth for this project because of forward application)
    
    properties
        
        val     % value
        der     % derivative
        
    end
    
    methods
        
        % Constructor
        % MUST be defined here, not in a separate file
        function obj = ADvar(a,b)
        
            % No input variables
            if nargin == 0
                
                obj.val = [];
                obj.der = [];
            
            % Only an input variable, considered constant
            elseif nargin == 1
                
                obj.val = a;
                obj.der = 0;

                % this may be useful but it requires testing
                % [r,c] = size(a);
                % obj.der = zeros(r, c);
            
            % Two input variables, considered as value and derivate
            else
                
                obj.val = a;
                obj.der = b;
                
            end
        end
        
        % Ad-hoc redefinition of some used functions
        % Only declarations here; definitions in separate files
        
        out = eq(u,v)
        
        out = ne(u,v)
        
        out = gt(u,v)
        
        out = ge(u,v)
        
        out = lt(u,v)
        
        out = le(u,v)
        
        u = uminus(u)
        
        u = cross(u,v)
        
        out = minus(u,v)
        
        out = plus(u,v)
        
        out = mtimes(u,v)
        
        out = times(u,v)
        
        out = mrdivide(u,v)
        
        out = rdivide(u,v)
        
        out = mpower(u,v)
        
        out = power(u,v)
        
        out = double(u)
        
        u = sqrt(u)
        
        u = sin(u)
        
        u = asin(u)
        
        u = cos(u)
        
        u = acos(u)
        
        u = tan(u)
        
        Y = atan2(Y,X)
        
        u = atan(u)
        
        u = log(u)
        
        u = exp(u)
        
        u = abs(u)
        
        out = getVal(u)
        
        out = getDer(u)
        
        out = horzcat(varargin)
        
        out = vertcat(varargin)
        
        A = transpose(A)
        
        u = subsref(u,s)
        
        u = subsasgn(u,s,val)
        
        u = sum(u, dim)
        
    end
end