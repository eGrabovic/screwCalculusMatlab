function Jb = bodyJacobian(ktObj, numNode)
% BODYJACOBIAN(ktObj, numNode) computes the body jacobian JB of 
%   NUMNODE w.r.t. the root in frame {numNode} (and thus its pole).
%   Because of how we computed the FWKin, the [B_0, j] matrix is 
%   just G_globalJ{i}.
%   We do not really compute the expTw(twist, q) here, instead we
%   borrow it from the FWKin computations.
%   The jacobian uses the FWKin joints values.

    subTree = ktObj.k{numNode};
    n = length(subTree);
    g = eye(4);
    Jb = zeros(6, n);
    
    for j = n : -1 : 1
        hl = ktObj.h(subTree(j)); % helicoidal lead to identify joint type

        if hl == 0 % revolute
            X = [0; 0; 0; 0; 0; 1];

        elseif hl == Inf % prismatic
            X = [0; 0; 1; 0; 0; 0];

        else
            X = [0; 0; hl; 0; 0; 1];
        end

        g = ktObj.G_localJ{subTree(j)} * g;
        Jb(:, j) = adjoint(hom_mat_inv(g)) * X; % adjoint(g_j,n^-1) * Y_i
        g = ktObj.Goff(:, :, subTree(j)) * g;
    
    end
end