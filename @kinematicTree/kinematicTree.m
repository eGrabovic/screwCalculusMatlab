classdef kinematicTree < handle
% Class used to build and manage a mechanical Tree related to a kinematic
%   object, i.e. a robot arm.
%
% ktObj = kinematicTree(p, s, k, mu, ni, gT, h)

    properties
        
        Nb         % Number of bodies
        p          % Predecessor array (parent nodes of i-th joint with 
                   %                    polarity)
        s          % Successors array (children nodes of i-th joint with 
                   %                    polarity)
        lambda     % Parent array
        k          % Cell array of joints from base to node; 
                   %    i -> element of the array
        mu         % Cell array of direct childrens of node i
        ni         % Cell array of all the nodes cointained in the sub-tree
                   %    with node i as root
        leaves     % Array of all the nodes without childrens
        Goff       % Array of offset homogenous matrices from joint i-1 to 
                   %    joint i (both with respect to body lambda(i) )
        h          % Array of helical lead of each joint 
                   %    (0 -> revolute; Inf-> prismatic; otherwise helicoidal)
        I          % Inertias
        G_local    % Local transform matrices for each body
        G_localJ   % Local transform matrices for each joint
        G_global   % Global transform matrices for each body
        G_globalJ  % Global transform matrices for each joint
        V          % Twists
        Vdot       % Twist derivative
        
        
        % ***!!!  mu{i} and ni{i} are associated to the (i-1)-th node !!!***
    end
    
    properties (Hidden = true)
        
        graphics
        adjoints
        
    end
    
    methods (Static, Hidden = true)
        
        [x, y, z] = createCylinder(r, len)
        
        points = createParallelepiped(edge, len)
        
    end
    
    methods (Static, Access = private)

        g = jointFun(x, h)
        
    end
    
    methods
        
        % Constructor
        % Must be defined here, not in another file
        function ktObj = kinematicTree(p, s, gT, h)
            
            ktObj.p = p;
            ktObj.s = s;
            ktObj.lambda = min(p, s);
            ktObj.Goff = gT;
            ktObj.h = h;
            ktObj.Nb = length(p);
            
            % compute the i-th body sub tree from node to body
            for i = 1:ktObj.Nb
                c = p(i);
                ktObj.k{i} = i;
                while c ~= 0 % until we arrive to root
                    ktObj.k{i} = [c, ktObj.k{i}];
                    if c == p(c) % or until the predecessor is the node itself
                        break
                    end
                    c = p(c);
                end
            end
            
            % compute the i-th body direct childrens
            for i = 1:ktObj.Nb+1
               ktObj.mu{i} = find(ktObj.lambda == i-1); 
            end
            
            for i = 1:ktObj.Nb+1 % compute i-th body sub tree with all the childerns
                index = i;
                ktObj.ni{i} = i - 1;
                while ~isempty(index)
                    ktObj.ni{i} = [ktObj.ni{i} ktObj.mu{index}];
                    index = [ktObj.mu{index}] + 1;
                end
            end
            
            ktObj.leaves = nan(1, ktObj.Nb + 1);
            cont = 0;
            for i = 1:length(ktObj.mu)
                if isempty(ktObj.mu{i})
                    cont = cont+1;
                    ktObj.leaves(cont) = i-1;
                end
            end
            
        end

        % Ad-hoc redefinition of some used functions
        % Only signatures here; definitions in separate files

        FWKin(ktObj, q)
        
        Jb = bodyJacobian(ktObj, numNode)
        
        plotInit(ktObj, varargin)
        
        updatePlot(ktObj)
        
        [tau, F] = RNEAdyn(ktObj, q, qd, qdd, V0, V0d, Fext)
        
        qdd = ABAdyn(ktObj, q, qd, tau, V0, V0d, Fext)
        
        qdd = ABAdynCasadi(ktObj, q, qd, tau, V0, V0d, Fext)
        
    end
    

end