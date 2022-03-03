function qdd = ABAdynCasadi(ktObj, q, qd, tau, V0, V0d, Fext)
% TODO
% q:    istantaneous joint configuration
% qd:   istantaneous joint velocity
% V0:   istantaneous base velocity
% V0d:  istantaneous base acceleration
% tau:  istantaneous joint action (forces/moments)
    
    cl = class(q);
    % Pose computation
    ktObj.FWKin(q);
    
    % Twist computation
    twists = nan(6, ktObj.Nb + 1, cl);
    twists(:, 1) = V0;
    X = nan(6, ktObj.Nb, cl);
    AdGinv = cell(1, ktObj.Nb);
    a      = cell(1, ktObj.Nb);

    %% Forward propagation of twists
    for i = 1:ktObj.Nb
        hl = ktObj.h(i); % helicoidal lead to identify joint type
        if hl == 0 % revolute
            X(:, i) = [0; 0; 0; 0; 0; 1];

        elseif hl == Inf % prismatic
            X(:, i) = [0; 0; 1; 0; 0; 0];

        else % helical
            X(:, i) = [0; 0; hl; 0; 0; 1];
        end

        AdGinv{i} = adjointInv(ktObj.G_local{i});
        twists(:, i+1) = AdGinv{i} * twists(:, ktObj.lambda(i) + 1) + ...
                         X(:, i) .* qd(i);
        a{i} = ad(X(:, i)) * qd(i);
    end

    %% Backward propagation of projected inertia and biases
    % Initialization of inertias
    Mtilde = cell(1, ktObj.Nb);
    Mbar = cell(1, ktObj.Nb);
    btilde = nan(6, ktObj.Nb, cl);
    bi     = nan(6, ktObj.Nb, cl);
    
    for i = ktObj.Nb:-1:1
        % Heavy math computations, see literature for documentation
        % (Meccanica dei Robot - Gabiccini ABA)
        Mtilde{i} = ktObj.I{i}; % Mtilde = M_i as initialization
        bi(:, i) = adStar(twists(:, i+1)) * ktObj.I{i} * twists(:, i+1);
        btilde(:, i) = bi(:, i) - Fext(:, i);
        m = ktObj.mu{i+1};

        for j = 1:length(m)
           children = m(j);
           ADstar = AdGinv{children}.';
           Mtilde{i} =  Mtilde{i} + ADstar * Mbar{children} * AdGinv{children};
           btilde(:, i) =  btilde(:, i) + ADstar * btilde(:, children) - ...
               ADstar * Mbar{children} * a{children} * AdGinv{children} * twists(:, i+1)+...
               (1 ./ (X(:, children).' * Mtilde{children} * X(:, children))) * ...
               ADstar * Mtilde{children} * X(:, children) * ...
               (tau(children) - X(:, children).' * btilde(:, children));
        end

        Mbar{i} = (eye(6) - Mtilde{i} * X(:,i) * X(:,i).' ./ ...
                  ((X(:, i)).' * Mtilde{i} * X(:,i))) * Mtilde{i};
    end
    
    %% Forward propagation of accelerations
    twistsD = nan(6, ktObj.Nb + 1, cl);
    twistsD(:,1) = V0d;
    qdd = nan(1, ktObj.Nb, cl);

    for i = 1:ktObj.Nb
        Aji = AdGinv{i};

        % i+1-th twist actually refers to the i-th body
        Vlambda = Aji * twists(:, ktObj.lambda(i) + 1);
        VlambdaD = Aji * twistsD(:, ktObj.lambda(i) + 1);
        qdd(i) = (tau(i) - X(:, i).' * (Mtilde{i} * (VlambdaD - a{i} * Vlambda) + ...
                 btilde(:, i))) ./ (X(:, i).' * Mtilde{i} * X(:, i));
        twistsD(:, i+1) = VlambdaD + X(:, i) .* qdd(i) - a{i} * Vlambda;
    end
end