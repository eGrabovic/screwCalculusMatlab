function FWKin(ktObj, q)
% FWKIN(ktObj, q) computes the forward kinematics with the local POE
%   formulation. All the intermediate transform matrices are stored.

    for i = 1:ktObj.Nb
        % twist exponential joint formulation
        ktObj.G_localJ{i} = kinematicTree.jointFun(q(i), ktObj.h(i));
        ktObj.G_local{i} = ktObj.Goff(:,:,i) * ktObj.G_localJ{i};

        if ktObj.lambda(i) == 0
            ktObj.G_global{i} = ktObj.G_local{i};
            ktObj.G_globalJ{i} = ktObj.Goff(:,:,i);
        else
            % property that lambda(i) < i
            ktObj.G_global{i} = ktObj.G_global{ktObj.lambda(i)} * ...
                                ktObj.G_local{i};
            ktObj.G_globalJ{i} = ktObj.G_global{ktObj.lambda(i)} * ...
                                 ktObj.Goff(:,:,i);
        end
    end
end