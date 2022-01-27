classdef kinematicTree < handle
%
% obj = kinematicTree(p, s, k, mu, ni, gT, h)
%
% TODO: help for every member function
    properties
        
        Nb         % number of bodies
        p          % predecessor array (parent nodes of i-th joint with polarity)
        s          % successors array (children nodes of i-th joint with polarity)
        lambda     % parent array
        k          % cell array of joints from base to node ; i -> element of the array
        mu         % cell array of direct childrens of node i
        ni         % cell array of all the nodes cointained in the sub-tree with node i as root
        leaves     % array of all the nodes without childrens
        Goff       % array of offset homogenous matrices from joint i-1 to joint i (both with respect to body lambda(i) )
        h          % array of helical lead of each joint (0 -> revolute; Inf-> prismatic; otherwise helical joint
        I          % inertias
        G_local    % local transform matrices for each body
        G_localJ   % local transform matrices for each joint
        G_global   % global transform matrices for each body
        G_globalJ  % global transform matrices for each joint
        V          % twists
        Vdot       % twist derivative
        
        
        % ***!!!  mu{i} and ni{i} are associated to the (i-1)-th node !!!***
    end
    
    properties (Hidden = true)
        
        graphics
        adjoints
        
    end
    
    methods (Static, Hidden = true)
        
        function [x, y, z] = createCylinder(r, len)
            
            theta = linspace(0, 2*pi, 15);
            x = r.*cos(theta);
            y = r.*sin(theta);
            z = len/2.*ones(1,15);
            
        end
        
        function points = createParallelepiped(edge, len)
            
            a = edge;
            if length(len) == 2
                h1 = len(1);
                h2 = len(2);
            else
                h1 = 0;
                h2 = len;
            end
            P1 = [a/2;-a/2;-h1]; P2 = [a/2;a/2;-h1]; P3 = [a/2;a/2;h2]; P4 = [a/2;-a/2;h2];
            P5 = [-a/2;a/2;-h1]; P6 = [-a/2;a/2;h2]; P7 = [-a/2;-a/2;-h1]; P8 = [-a/2;-a/2;h2];
            points = [P1,P2,P3,P4,P5,P6,P7,P8];
            
        end
        
    end
    
    methods (Static, Access = private)
        
        function g = jointFun(x, h)
            
            switch h
                case 0
                    g = [cos(x) -sin(x) 0    0;
                         sin(x)  cos(x) 0    0;
                         0       0      1    0;
                         0       0      0    1];
                case Inf
                    g = [eye(3), [0; 0; x]; 0 0 0 1];
                otherwise
                    g = [cos(x) -sin(x)  0     0;
                         sin(x)  cos(x)  0     0;
                          0       0      1  x.*h;
                          0       0      0     1];
            end
            
        end
    end
    
    methods
        
        function obj = kinematicTree(p, s, gT, h)
            
            obj.p = p;
            obj.s = s;
            obj.lambda = min(p, s);
            obj.Goff = gT;
            obj.h = h;
            obj.Nb = length(p);
            
            % compute the i-th body sub tree from node to body
            for i = 1:obj.Nb
                c = p(i);
                obj.k{i} = i;
                while c ~= 0 % until we arrive to root
                    obj.k{i} = [c, obj.k{i}];
                    if c == p(c) % or until the predecessor is the node itself
                        break
                    end
                    c = p(c);
                end
            end
            
            % compute the i-th body direct childrens
            for i = 1:obj.Nb+1
               obj.mu{i} = find(obj.lambda == i-1); 
            end
            
            for i = 1:obj.Nb+1 % compute i-th body sub tree with all the childerns
                index = i;
                obj.ni{i} = i-1;
                while ~isempty(index)
                    obj.ni{i} = [obj.ni{i} obj.mu{index}];
                    index = [obj.mu{index}] + 1;
                end
            end
            
            obj.leaves = nan(1, obj.Nb+1);
            cont = 0;
            for i = 1:length(obj.mu)
                if isempty(obj.mu{i})
                    cont = cont+1;
                    obj.leaves(cont) = i-1;
                end
            end
            
        end
        
        function FWKin(this, q)
        % FWKIN(this, q) computes the forward kinematics with the local POE
        %   formulation. All the intermediate transform matrices are stored.

            for i = 1:this.Nb
             this.G_localJ{i} = kinematicTree.jointFun(q(i), this.h(i)); % twist exponential joint formulation
             this.G_local{i} = this.Goff(:,:,i)*this.G_localJ{i};
             if this.lambda(i) == 0
                this.G_global{i} = this.G_local{i};
                this.G_globalJ{i} = this.Goff(:,:,i);
             else
                this.G_global{i} = this.G_global{this.lambda(i)}*this.G_local{i}; % property that lambda(i) < i
                this.G_globalJ{i} = this.G_global{this.lambda(i)}*this.Goff(:,:,i);
             end
            end
        end
        
        function Jb = bodyJacobian(this, numNode)
        % BODYJACOBIAN(this, numNode) computes the body jacobian JB of 
        %   NUMNODE w.r.t. the root in frame {numNode} (and thus its pole).
        %   Because of how we computed the FWKin, the [B_0, j] matrix is 
        %   just G_globalJ{i}.
        %   We do not really compute the expTw(twist, q) here, instead we
        %   borrow from the FWKin computations.
        %   The jacobian uses the FWKin joints values.
            
            subTree = this.k{numNode};
            n = length(subTree);
            g = eye(4);
            Jb = zeros(6, n);
            
            for j = n : -1 : 1
                hl = this.h(subTree(j)); % helical lead to identify joint type
                if hl == 0 % revolute
                    X = [0;0;0;0;0;1];
                elseif hl == Inf % prismatic
                    X = [0;0;1;0;0;0];
                else
                    X = [0;0;hl;0;0;1];
                end
                g = this.G_localJ{subTree(j)}*g;
                Jb(:, j) = adjoint(hom_mat_inv(g))*X; % adjoint(g_j,n^-1)*Y_i
                g = this.Goff(:,:,subTree(j))*g;
                
            end
        end
        
%         function spatialJacobian(this, numNode)
%             % spatial jacobion of numNode w.r.t. to the root in root's
%             % frame
%             
%             
%         end
        
        function plotInit(this, varargin)
        %
        %
        %
            
            if isempty(varargin)
                figure('color', 'w'); hold on; axis equal
                this.graphics.axes = gca;
            else
                this.graphics.axes = varargin{1};
            end
            
            % precomputing end effector graphics points
            S(:,1) = [0;0;0];
            S(:,2) = [0;0;0.5];
            S(:,3) = [0.5;0;0.5];
            S(:,4) = [0.5;0;0.7];
            
            % precomputing revolute joint gaphics
            [x, y, z] = kinematicTree.createCylinder(0.1, 0.5);
            ptsR = [x;y;z];
            % precomputing prismatic joint gaphics
            ptsP = kinematicTree.createParallelepiped(0.1, [0.25 0.25]);
            ptsPslide = kinematicTree.createParallelepiped(0.05, [0.4 0.4]);
            
            for i = 1:this.Nb % joint plots
                % initi the hgtr handles
                this.graphics.homogTransfJoint{i} = hgtransform(Parent = this.graphics.axes);
                this.graphics.homogTransfBody{i} = hgtransform(Parent = this.graphics.axes);
                this.graphics.homogTransfJointNext{i} = hgtransform(Parent = this.graphics.axes);
                
                % update the hgtr handles
                this.graphics.homogTransfJoint{i}.Matrix = this.G_globalJ{i};
                this.graphics.homogTransfJointNext{i}.Matrix = this.G_global{i};
                if this.lambda(i) == 0
                    this.graphics.homogTransfBody{i}.Matrix = eye(4);
                else
                    this.graphics.homogTransfBody{i}.Matrix = this.G_global{this.lambda(i)};
                end
                
                if this.h(i) ~= Inf % revolute joint or helical
                    
                    surf([ptsR(1,:); ptsR(1,:)], [ptsR(2,:); ptsR(2,:)], [ptsR(3,:); -ptsR(3,:)], 'facecolor', 'r', 'edgecolor', 'none', 'Parent', this.graphics.homogTransfJoint{i});
                    fill3(ptsR(1,:), ptsR(2,:), ptsR(3,:), 'r', 'Parent', this.graphics.homogTransfJoint{i});
                    fill3(ptsR(1,:), ptsR(2,:), -ptsR(3,:), 'r', 'Parent', this.graphics.homogTransfJoint{i});
                else   % prismatic joint 
                    
                    faces = [1,2,3,4;2,5,6,3; 3,6,8,4; 8,6,5,7; 1,4,8,7; 1,7,5,2];
                    patch('faces', faces, 'vertices', ptsP.',  'facecolor', 'green', 'Parent', this.graphics.homogTransfJoint{i});
                    patch('faces', faces, 'vertices', ptsPslide.',  'facecolor', 'green', 'Parent', this.graphics.homogTransfJointNext{i});
                end
                
                P = this.Goff(:,:,i);
                P = P(1:3, 4);
                L(:, 1) = [0;0;0];
                L(:, 2) = [P(1);0;0];
                L(:, 3) = [P(1);P(2);0];
                L(:, 4) = [P(1);P(2);P(3)];
                line(L(1,:), L(2,:), L(3,:), 'color', 'k', 'linewidth', 1.4, 'Parent', this.graphics.homogTransfBody{i});
                
                if isempty(this.mu{i+1}) % plot also a placeholder end effector if the i-th body has no children
                   line(S(1,:), S(2,:), S(3,:), 'color', 'b', 'linewidth', 1.4, 'Parent', this.graphics.homogTransfJointNext{i});
                end
            end
            
        end
        
        function updatePlot(this)
        %
        % 
        %
            
            % check if plot initialization still exists
            if isempty(this.graphics) || ~isvalid(this.graphics.axes)
                error('no current plot to update exists, initialize the plot first with "treeObject.plotInit" ');
            end
            
            % updating the kinematics
%             FWKin(this, q);
            
            for i = 1:this.Nb % joint plots
                % just update the hgtransform handles and we're done
                this.graphics.homogTransfJoint{i}.Matrix = this.G_globalJ{i};
                this.graphics.homogTransfJointNext{i}.Matrix = this.G_global{i};
                if this.lambda(i) == 0
                    this.graphics.homogTransfBody{i}.Matrix = eye(4);
                else
                    this.graphics.homogTransfBody{i}.Matrix = this.G_global{this.lambda(i)};
                end
                
            end
        end
        
        function [tau, F] = RNEAdyn(this, q, qd, qdd, V0, V0d, Fext)
        %
        %
        %
            
            %% forward recursive computation
            % posture computation
            this.FWKin(q);
            
            % twist computation
            twists = nan(6, this.Nb+1);
            twists(:, 1) = V0;
            twistsd(:, 1) = V0d;
            X = nan(6, this.Nb);
            
            for i = 1:this.Nb
                hl = this.h(i); % helical lead to identify joint type
                if hl == 0 % revolute
                    X(:,i) = [0;0;0;0;0;1];
                elseif hl == Inf % prismatic
                    X(:,i) = [0;0;1;0;0;0];
                else % helical
                    X(:,i) = [0;0;hl;0;0;1];
                end
                ADjG = adjointInv(this.G_local{i});
                ADjGV = ADjG*twists(:,this.lambda(i) + 1);
                Xd = X(:,i).*qd(i);
                twists(:, i+1) = ADjGV + Xd;
                twistsd(:, i+1) = X(:,i).*qdd(i) + ADjG*twistsd(:,this.lambda(i) + 1) + ad(ADjGV)*Xd;
            end
            %% backward recursive computation
            
            % generalized forces initialization
            F = nan(6, this.Nb);
            % actuator component
            tau = nan(1, this.Nb);
 
            for i = this.Nb:-1:1
                % computing the forces on all the direct childrens
                Fch = zeros(6,1);
                m = this.mu{i+1};
                if ~isempty(m)
                    for j = 1:length(m)
                        Fch = Fch + adjointStar(this.G_local{m(j)})*F(:, m(j));
                    end
                end
                II = this.I{i}*twistsd(:, i+1) + adStar(twists(:,i+1))*this.I{i}*twists(:, i+1); % inertial force
                F(:, i) = Fch - Fext(:, i) + II;  % equilibrium
                tau(i) = X(:, i).'*F(:, i);       % extracting the joint component
            end
        end
        
        function qdd = ABAdyn(this, q, qd, tau, V0, V0d, Fext)
        %
        % q:    istantaneous joint configuration
        % qd:   istantaneous joint velocity
        % V0:   istantaneous base velocity
        % V0d:  istantaneous base acceleration
        % tau:  istantaneous joint action (forces/moments)
            
%             cl = class(q);
            % compute posture
            this.FWKin(q);
            
            % twist computation
            twists = nan(6, this.Nb+1);
            twists(:, 1) = V0;
            X = nan(6, this.Nb);
            AdGinv = nan(6,6,this.Nb);
            a      = nan(6, 6, this.Nb);
            % forward propagation of twists
            for i = 1:this.Nb
                hl = this.h(i); % helical lead to identify joint type
                if hl == 0 % revolute
                    X(:,i) = [0;0;0;0;0;1];
                elseif hl == Inf % prismatic
                    X(:,i) = [0;0;1;0;0;0];
                else % helical
                    X(:,i) = [0;0;hl;0;0;1];
                end
                AdGinv(:,:,i) = adjointInv(this.G_local{i});
                twists(:, i+1) = AdGinv(:,:,i)*twists(:,this.lambda(i) + 1) + X(:,i).*qd(i);
                a(:,:,i) = ad(X(:,i))*qd(i);
            end

            % backward propagation of projected inertia and biases
            
            % intialization of inertias
            Mtilde = nan(6, 6, this.Nb);
            Mbar   = nan(6, 6, this.Nb);
            btilde = nan(6, this.Nb);
            bi     = nan(6, this.Nb);
            
            for i = this.Nb:-1:1
                % heavy math computations, see literature for documentation
                % (meccanica dei robot Gabiccini ABA)
                Mtilde(:,:,i) = this.I{i}; % Mtilde = M_i as initialization
                bi(:,i) = adStar(twists(:,i+1))*this.I{i}*twists(:,i+1);
                btilde(:,i) = bi(:,i) - Fext(:, i);
                m = this.mu{i+1};
                for j = 1:length(m)
                   children = m(j);
                   ADstar = AdGinv(:,:,children).';
                   Mtilde(:,:,i) =  Mtilde(:,:,i) + ADstar*Mbar(:,:,children)*AdGinv(:,:,children);
                   btilde(:,i) =  btilde(:,i) + ADstar*btilde(:,children) - ADstar*Mbar(:,:,children)*a(:, :, children)*AdGinv(:,:,children)*twists(:, i+1)+...
                       (1./(X(:,children).'*Mtilde(:,:,children)*X(:, children)))*ADstar*Mtilde(:,:,children)*X(:,children)*(tau(children) - X(:,children).'*btilde(:, children));
                end
                Mbar(:,:,i) = (eye(6) - Mtilde(:,:,i)*X(:,i)*X(:,i).'./((X(:,i)).'*Mtilde(:,:,i)*X(:,i)))*Mtilde(:,:,i);
            end
            
            % forward propagation of accelerations
            twistsD = nan(6, this.Nb+1);
            twistsD(:,1) = V0d;
            qdd = nan(1,this.Nb);
            for i = 1:this.Nb
                Aji = AdGinv(:,:,i);
                Vlambda = Aji*twists(:,this.lambda(i)+1); % i+1-th twist actually refers to the i-th body
                VlambdaD = Aji*twistsD(:,this.lambda(i)+1);
                qdd(i) = (tau(i) - X(:,i).'*(Mtilde(:,:,i)*(VlambdaD - a(:,:,i)*Vlambda) + btilde(:,i)))./(X(:,i).'*Mtilde(:,:,i)*X(:,i));
                twistsD(:, i+1) = VlambdaD + X(:,i).*qdd(i) - a(:,:,i)*Vlambda;
            end
        end
        
        function qdd = ABAdynCasadi(this, q, qd, tau, V0, V0d, Fext)
        %
        % q:    istantaneous joint configuration
        % qd:   istantaneous joint velocity
        % V0:   istantaneous base velocity
        % V0d:  istantaneous base acceleration
        % tau:  istantaneous joint action (forces/moments)
            
            cl = class(q);
            % compute posture
            this.FWKin(q);
            
            % twist computation
            twists = nan(6, this.Nb+1, cl);
            twists(:, 1) = V0;
            X = nan(6, this.Nb, cl);
            AdGinv = cell(1, this.Nb);
            a      = cell(1, this.Nb);
            % forward propagation of twists
            for i = 1:this.Nb
                hl = this.h(i); % helical lead to identify joint type
                if hl == 0 % revolute
                    X(:,i) = [0;0;0;0;0;1];
                elseif hl == Inf % prismatic
                    X(:,i) = [0;0;1;0;0;0];
                else % helical
                    X(:,i) = [0;0;hl;0;0;1];
                end
                AdGinv{i} = adjointInv(this.G_local{i});
                twists(:, i+1) = AdGinv{i}*twists(:,this.lambda(i) + 1) + X(:,i).*qd(i);
                a{i} = ad(X(:,i))*qd(i);
            end

            % backward propagation of projected inertia and biases
            
            % intialization of inertias
            Mtilde = cell(1, this.Nb);
            Mbar = cell(1, this.Nb);
            btilde = nan(6, this.Nb, cl);
            bi     = nan(6, this.Nb, cl);
            
            for i = this.Nb:-1:1
                % heavy math computations, see literature for documentation
                % (meccanica dei robot Gabiccini ABA)
                Mtilde{i} = this.I{i}; % Mtilde = M_i as initialization
                bi(:,i) = adStar(twists(:,i+1))*this.I{i}*twists(:,i+1);
                btilde(:,i) = bi(:,i) - Fext(:, i);
                m = this.mu{i+1};
                for j = 1:length(m)
                   children = m(j);
                   ADstar = AdGinv{children}.';
                   Mtilde{i} =  Mtilde{i} + ADstar*Mbar{children}*AdGinv{children};
                   btilde(:,i) =  btilde(:,i) + ADstar*btilde(:,children) - ADstar*Mbar{children}*a{children}*AdGinv{children}*twists(:, i+1)+...
                       (1./(X(:,children).'*Mtilde{children}*X(:, children)))*ADstar*Mtilde{children}*X(:,children)*(tau(children) - X(:,children).'*btilde(:, children));
                end
                Mbar{i} = (eye(6) - Mtilde{i}*X(:,i)*X(:,i).'./((X(:,i)).'*Mtilde{i}*X(:,i)))*Mtilde{i};
            end
            
            % forward propagation of accelerations
            twistsD = nan(6, this.Nb+1, cl);
            twistsD(:,1) = V0d;
            qdd = nan(1,this.Nb, cl);
            for i = 1:this.Nb
                Aji = AdGinv{i};
                Vlambda = Aji*twists(:,this.lambda(i)+1); % i+1-th twist actually refers to the i-th body
                VlambdaD = Aji*twistsD(:,this.lambda(i)+1);
                qdd(i) = (tau(i) - X(:,i).'*(Mtilde{i}*(VlambdaD - a{i}*Vlambda) + btilde(:,i)))./(X(:,i).'*Mtilde{i}*X(:,i));
                twistsD(:, i+1) = VlambdaD + X(:,i).*qdd(i) - a{i}*Vlambda;
            end
        end
        
    end
    

end