function [tau, F] = RNEAdyn(ktObj, q, qd, qdd, V0, V0d, Fext)
% TODO
%
%
    
    %% Forward recursive computation
    % Pose computation
    ktObj.FWKin(q);
    
    % Twist computation
    twists = nan(6, ktObj.Nb + 1);
    twists(:, 1) = V0;
    twistsd(:, 1) = V0d;
    X = nan(6, ktObj.Nb);
    
    for i = 1:ktObj.Nb
        hl = ktObj.h(i); % helicoidal lead to identify joint type

        if hl == 0 % revolute
            X(:,i) = [0; 0; 0; 0; 0; 1];

        elseif hl == Inf % prismatic
            X(:,i) = [0; 0; 1; 0; 0; 0];

        else % helicoidal
            X(:,i) = [0; 0; hl; 0; 0; 1];
        end

        ADjG = adjointInv(ktObj.G_local{i});
        ADjGV = ADjG*twists(:,ktObj.lambda(i) + 1);
        Xd = X(:,i) .* qd(i);
        twists(:, i+1) = ADjGV + Xd;
        twistsd(:, i+1) = X(:,i) .* qdd(i) + ...
                          ADjG * twistsd(:,ktObj.lambda(i) + 1) + ...
                          ad(ADjGV) * Xd;
    end
    %% Backward recursive computation
    % Generalized forces initialization
    F = nan(6, ktObj.Nb);

    % Actuator component
    tau = nan(1, ktObj.Nb);

    for i = ktObj.Nb:-1:1
        % Computing the forces on all the direct childrens
        Fch = zeros(6, 1);
        m = ktObj.mu{i+1};

        if ~isempty(m)
            for j = 1:length(m)
                Fch = Fch + adjointStar(ktObj.G_local{m(j)}) * F(:, m(j));
            end
        end

        % Inertial force
        II = ktObj.I{i} * twistsd(:, i+1) + ...
             adStar(twists(:,i+1)) * ktObj.I{i} * twists(:, i+1); 
        % Equilibrium
        F(:, i) = Fch - Fext(:, i) + II;
        % Extracting the joint component
        tau(i) = X(:, i).' * F(:, i);
    end
end