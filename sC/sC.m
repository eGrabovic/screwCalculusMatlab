classdef sC
% Class that serves as "library" (pseudo screwCalculus) for different
% functions useful to compute the kinematics of a shaping machine for
% spiral bevel wheels.
%
%       FUNCTIONS (callable as METHODS: sC.method() )
%
%   N.B. You can improve calling speed of library functions renaming the
%   class in the workspace (i.e. SC = sC; SC.hat() calls the hat() function)
%
% - hat(vector)
%
% - rotX(alpha) (the MATLAB one does not allow for symbolic computations)
%
% - rotY(alpha)
%
% - rotZ(alpha)
%
% - hom_mat(rotationMatrix, displacement)
%
% - hom_mat_inv(Gst): inverse of a homog. matrix transf. Gst
%
% - adjoint(rotationMatrix,displacement) or adjoint(homogeneousmatrix)
%
% - expSkew(axis,th) -> exp(omegahat*theta)
%
% - unitTwist(jointType,axis,q) (NO HELIC. JOINTS): reference unit Twist
%
% - expTw(unitTwist,th) ( NO HELIC. JOINS): twist Exponential
%
% - FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): forward Kinematics
%
% - AD_FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): automatic diff.
%   and Forward Kinematics (automatic diff. variable on joint vars).
%
% - SpJac({Y1,var1},{Y2,var2},...,{Yn,varn}): spatial Jacobian
%
% - BodyJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): body Jacobian
%
%   TODO: copy external help function or what?

    methods(Static)
        
        function X = hat(v)
        % HAT(vec) transforms a R3 column vector VEC in its 'hat' R3x3
        % antisymmetric matrix form, or transforms a R3 column twist vector in a
        % homogeneous 'hat' R4x4 matrix.  
        
            if all(size(v) ~= [3,1]) || all(size(v) ~= [6,1])
                error('v requires to be a 3 x 1 matrix (vector)')
            end
            
            if all(size(v) == [3,1])
                X = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
            end
            
            if all(size(v) == [6,1])
                vel = [v(1);v(2);v(3)];
                omeg = [v(4);v(5);v(6)];
                X = [sC.hat(omeg) vel;0 0 0 0];
            end
        end
        
        function [n, theta] = rotToAxisAngle(R)
            
           theta = acos((trace(R)-1)/2) ;
           n = 1/(2*sin(theta)).*[R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
            
        end
        
        function R = rotX(alpha)
            
            R = [1 0 0;0 cos(alpha) -sin(alpha);0 sin(alpha) cos(alpha)];
            
        end
        
        function R = rotY(alpha)
            %
            % funzione matrice di rotazione attorno asse Y per acconsentire
            % il calcolo simbolico
            %
            R = [cos(alpha) 0 sin(alpha);0 1 0; -sin(alpha) 0 cos(alpha)];
        end
        
        function R = rotZ(alpha)
            %
            % funzione matrice di rotazione attorno asse Z per acconsentire
            % il calcolo simbolico
            %
            R = [cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1];   
        end
        
        function R = rotZ2D(alpha)
            
            R = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
        end
        
        function R = diffRotZ2D(alpha)
        %
        %
        %
            R = [-sin(alpha) -cos(alpha);cos(alpha) -sin(alpha)];
           
        end
        
        function Gst = hom_mat(rotm,disp)
        % HOM_MAT(rotm, disp) returns the homogeneous rototranslational matrix OUT
        %   related to the rotation matrix ROTM and the displacement DISP.

            if all(size(rotm) ~= [3,3]) || all(size(disp) ~= [3,1])
                error('rotation matrix requires to have 3x3 dimension and the displacement vector requires to be a 3x1 matrix')
            end
            Gst = [rotm disp; 0 0 0 1];
        end
        
        function GstInv = hom_mat_inv(Gst)
        % HOM_MAT_INV(Gst) computes the inverse GSTINV of a homogeneous matrix GST.

            if all(size(Gst) ~= [4,4])
                error('Homogenous transformation matrix requires to have 4x4 dimension')
            end
            Rst = Gst(1:3,1:3);
            dst = Gst(1:3,4);
            GstInv = [Rst',-Rst'*dst;0,0,0,1];
        end
        
        function GstInv = hom_mat_inv2D(Gst)
        %
        %

            Rst = Gst(1:2,1:2);
            dst = Gst(1:2,3);
            GstInv = [Rst.',-Rst.'*dst;0,0,1];
            
            
        end
        
        function X = adjoint(rotm,disp)
            %
            % adjoint(Gst)
            % trasformata aggiunta di Gst. Matrice 6x6 che può trasformare
            % twist congruenti in base alla trasformazione omogenea di Gst.
            % adjoint(rotm,disp)
            % nel caso di 2 input vengono assegnate la matrice di rotazione
            % (rotm) e vett traslazione (disp) di Gst.
            %
            if nargin == 1
                if all(size(rotm) ~= [4,4])
                    error('size of Gst homogeneous matrix requires to be 4x4')
                end
                R = rotm(1:3,1:3);
                d = rotm(1:3,4);
                X = [R,sC.hat(d)*R;zeros(3,3),R];
            else
                if all(size(rotm) ~= [3,3]) || all(size(disp) ~= [3,1])
                    error('rotation matrix requires to be a 3x3 matrix and the displacement vector requires to be a 3x1 matrix')
                end
                R = rotm;
                d = disp;
                X = [R sC.hat(d)*R;zeros(3,3) R];
            end
        end
        
        function X = expSkew(axis,th)
            %
            % Calcola l'esponenziale e^(omega*theta);
            % omega = versore asse (axis);
            % theta = angolo (th);
            %
            % utilizza la formula di Rodriguez quando viene assegnato un
            % asse generico;
            %
            % vengono direttamente riportate invece le rotazioni elementari
            % nel caso l'asse coincida con uno degli assi principali (per
            % velocizzare i calcoli)
            if all(isnumeric(axis)) == true
                if abs(axis.'*axis) -1 >= 1e-8
                    error('axis requires to be a versor');
                end
                if abs(axis(1)) <= 1e-8 && abs(axis(2)) <= 1e-8
                    
                    if axis(3) > 0
                        
                        X = sC.rotZ(th);
                        
                    elseif axis(3) < 0
                        
                        X = sC.rotZ(-th);
                    end
                    
                elseif abs(axis(2)) <= 1e-8 && abs(axis(3)) <= 1e-8
                    
                    if axis(1) > 0
                        
                        X = sC.rotX(th);
                        
                    elseif axis(1) < 0
                        
                        X = sC.rotX(-th);
                    end
                    
                elseif abs(axis(1)) <= 1e-8 && abs(axis(3)) <= 1e-8
                    
                    if axis(2) > 0
                        
                        X = sC.rotY(th);
                        
                    elseif axis(2) < 0
                        
                        X = sC.rotY(-th);
                    end
                else
                    X = eye(3) + sC.hat(axis)*sin(th) + sC.hat(axis)*sC.hat(axis)*(1 - cos(th));
                    %rodriguez formula
                end
                
            else
                if abs(double(axis(1))) <= 1e-8 && abs(double(axis(2))) <= 1e-8
                    
                    if double(axis(3)) > 0
                        X = sC.rotZ(th);
                    elseif double(axis(3)) < 0
                        X = sC.rotZ(-th);
                    end
                elseif abs(double(axis(1))) <= 1e-8 && abs(double(axis(3))) <= 1e-8
                    
                    if double(axis(2)) > 0
                        X = sC.rotY(th);
                    elseif double(axis(2)) < 0
                        X = sC.rotY(-th);
                    end
                elseif abs(double(axis(2))) <= 1e-8 && abs(double(axis(3))) <= 1e-8
                    
                    if double(axis(1)) > 0
                        X = sC.rotX(th);
                    elseif double(axis(1)) < 0
                        X = sC.rotX(-th);
                    end
                else
                    X = eye(3) + sC.hat(axis)*sin(th) + sC.hat(axis)*sC.hat(axis)*(1 - cos(th));
                    % Rodriguez formula
                end
            end
        end
        
        function uT = unitTwist(jointType,axis,q)
            %
            % Crea un Twist unitario a partire dal tipo di giunto ('P' per
            % prismatico o 'R' per rotoidale), l'asse del giunto ed
            % eventualmente la distanza dall origine q se di tipo 'R'
            %
            if strcmpi(jointType,'P') == true
                uT = [axis;0;0;0];
            end
            
            if strcmpi(jointType,'R') == true
                uT = [-cross(axis,q);axis];
            end
            
        end
        
        function E = expTw(unitTwist,th)
            %
            % Calcolo dell esponenziale di un twist
            %
            % unitTwist = input twist unitario relativo al giunto
            %
            % th = variabile di giunto
            %
            if all(isnumeric(unitTwist)) == true
                if abs(unitTwist(4)) <= 1e-6 && abs(unitTwist(5)) <= 1e-6 && abs(unitTwist(6)) <= 1e-6
                    
                    axis = [unitTwist(1);unitTwist(2);unitTwist(3)];
                    R = eye(3);
                    d = axis*th;
                    
                else
                    
                    axis = [unitTwist(4);unitTwist(5);unitTwist(6)];
                    v = [unitTwist(1);unitTwist(2);unitTwist(3)];
                    q = -cross(v,axis);
                    R = sC.expSkew(axis,th);
                    d =(eye(3) - R)*q;
                    
                end
                E = sC.hom_mat(R,d);

            end
        end
        
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
            if n == 1
                gst = sC.expTw(varargin{1}{1},varargin{1}{2});
                gst = gst*gst0;
                return
            end
            gst = sC.expTw(varargin{1}{1},varargin{1}{2});
            for i = 2 : 1 : n
                gst = gst*sC.expTw(varargin{i}{1},varargin{i}{2});
            end
            gst = gst*gst0;
        end
        
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
                gst = ADexpTw(varargin{1}{1},varargin{1}{2});
                gst = gst*gst0;
                return
            end
            gst = ADexpTw(varargin{1}{1},varargin{1}{2});
            for i = 2 : 1 : n
                gst = gst*ADexpTw(varargin{i}{1},varargin{i}{2});
            end
            gst = gst*gst0;
        end
        
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
        
        function Jb = BodyJac(gst0,varargin)
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
            g = gst0;
            
            if n == 1
                Ytilde = sC.adjoint(sC.hom_mat_inv(g))*varargin{1}{1};
                Jb = Ytilde;
                return
            end
            
            Jb = zeros(6,n);
            
            for i = n : -1 : 1
                Ytilde = sC.adjoint(sC.hom_mat_inv(g))*varargin{i}{1};
                Jb(:,i) = Ytilde;
                g = sC.expTw(varargin{i}{1},varargin{i}{2})*g;
            end
            
        end
        
    end
end