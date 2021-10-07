function FKmat = ForwardKinematics(varargin) % DHFKine
    if (nargin < 1) || (nargin >2)
        return
    end
    
    % Initial homogeneous matrix
    H = eye(4);
    DHTableFull = varargin{1};
    DHTable = DHTableFull(:, 1:4); %str2double(DHTableFull(:,1:4));
    types = arrayfun(@string, DHTableFull(:,end)); % cell2mat(DHTableFull(:,end));
    n = size(DHTable, 1);

    if (nargin == 1)
        index = n;
    else
        index = varargin{2};
    end
    
    % Build Jacobian joint by joint
    for i = [1:index]
       H = H * FKRow(DHTable(i,:));
    end
       
    FKmat = H;
end