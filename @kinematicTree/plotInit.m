function plotInit(ktObj, varargin)
% TODO
%
%
    
    if isempty(varargin)
        figure('color', 'w'); hold on; axis equal
        ktObj.graphics.axes = gca;
    else
        ktObj.graphics.axes = varargin{1};
    end
    
    % Precomputing end-effector graphic points
    S(:,1) = [0; 0; 0];
    S(:,2) = [0; 0; 0.5];
    S(:,3) = [0.5; 0; 0.5];
    S(:,4) = [0.5; 0; 0.7];
    
    % Precomputing revolute joint gaphics
    [x, y, z] = kinematicTree.createCylinder(0.1, 0.5);
    ptsR = [x; y; z];

    % Precomputing prismatic joint gaphics
    ptsP = kinematicTree.createParallelepiped(0.1, [0.25 0.25]);
    ptsPslide = kinematicTree.createParallelepiped(0.05, [0.4 0.4]);
    
    % Joints plot
    for i = 1:ktObj.Nb
        % Initialize the hgtr handles
        ktObj.graphics.homogTransfJoint{i} = hgtransform(Parent = ktObj.graphics.axes);
        ktObj.graphics.homogTransfBody{i} = hgtransform(Parent = ktObj.graphics.axes);
        ktObj.graphics.homogTransfJointNext{i} = hgtransform(Parent = ktObj.graphics.axes);
        
        % Update the hgtr handles
        ktObj.graphics.homogTransfJoint{i}.Matrix = ktObj.G_globalJ{i};
        ktObj.graphics.homogTransfJointNext{i}.Matrix = ktObj.G_global{i};

        if ktObj.lambda(i) == 0
            ktObj.graphics.homogTransfBody{i}.Matrix = eye(4);
        else
            ktObj.graphics.homogTransfBody{i}.Matrix = ktObj.G_global{ktObj.lambda(i)};
        end
        
        if ktObj.h(i) ~= Inf % revolute or helicoidal joint
            surf([ptsR(1,:); ptsR(1,:)], [ptsR(2,:); ptsR(2,:)], [ptsR(3,:); -ptsR(3,:)], ...
                'facecolor', 'r', 'edgecolor', 'none', ...
                'Parent', ktObj.graphics.homogTransfJoint{i});
            fill3(ptsR(1,:), ptsR(2,:), ptsR(3,:), 'r', ...
                'Parent', ktObj.graphics.homogTransfJoint{i});
            fill3(ptsR(1,:), ptsR(2,:), -ptsR(3,:), 'r', ...
                'Parent', ktObj.graphics.homogTransfJoint{i});

        else % prismatic joint
            faces = [1,2,3,4;2,5,6,3; 3,6,8,4; 8,6,5,7; 1,4,8,7; 1,7,5,2];
            patch('faces', faces, 'vertices', ptsP.', 'facecolor', 'green',...
                  'Parent', ktObj.graphics.homogTransfJoint{i});
            patch('faces', faces, 'vertices', ptsPslide.', ...
                  'facecolor', 'green', ...
                  'Parent', ktObj.graphics.homogTransfJointNext{i});
        end
        
        P = ktObj.Goff(:,:,i);
        P = P(1:3, 4);
        L(:, 1) = [0; 0; 0];
        L(:, 2) = [P(1); 0; 0];
        L(:, 3) = [P(1); P(2); 0];
        L(:, 4) = [P(1); P(2); P(3)];
        line(L(1, :), L(2, :), L(3, :), 'color', 'k', 'linewidth', 1.4, ...
             'Parent', ktObj.graphics.homogTransfBody{i});
        
        % Plot also a placeholder end-effector if the i-th body has no children
        if isempty(ktObj.mu{i+1}) 
           line(S(1,:), S(2,:), S(3,:), 'color', 'b', 'linewidth', 1.4, ...
               'Parent', ktObj.graphics.homogTransfJointNext{i});
        end
    end
    
end