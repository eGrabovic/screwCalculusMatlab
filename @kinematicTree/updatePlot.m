function updatePlot(ktObj)
% TODO
% 
%
    
    % Check if plot initialization still exists
    if isempty(ktObj.graphics) || ~isvalid(ktObj.graphics.axes)
        error(['No current plot to update exists, initialize the plot ' ...
            'with "treeObject.plotInit" before using this function']);
    end
    
    % Update the kinematics
%             FWKin(ktObj, q);
    
    % Joints plot
    for i = 1:ktObj.Nb
        % Just update the hgtransform handles and we're done
        ktObj.graphics.homogTransfJoint{i}.Matrix = ktObj.G_globalJ{i};
        ktObj.graphics.homogTransfJointNext{i}.Matrix = ktObj.G_global{i};
        
        if ktObj.lambda(i) == 0
            ktObj.graphics.homogTransfBody{i}.Matrix = eye(4);
        else
            ktObj.graphics.homogTransfBody{i}.Matrix = ktObj.G_global{ktObj.lambda(i)};
        end
        
    end
end