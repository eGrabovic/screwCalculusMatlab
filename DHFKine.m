function FKmat = DHFKine(varargin)
% DHFKINE builds the Forward Kinematics Homogeneous matrix.
%
%     The main argument to pass is the Denavit-Hartemberg table, but the 
%     function handles up to 4 arguments, adapting to the situation:
%     1) only the table is supplied -> build entire Forward Kinematics
%     2) table and an index, specifying the joint number at which to stop
%     3) table, pre- and post- moltiplication matrix, for additional
%        transformations after the (whole) Forward Kinematics
%     4) all of the above


    % Check number of arguments
    if (nargin < 1) || (nargin >4)
        return
    end
    
    % Initialize and check first argument
    DHTableFull = varargin{1};
    DHTable = DHTableFull(:, 1:4);
    n = size(DHTable, 1);

    % Handle numeric DHTable converted in string because of the 'type'
    % column; non-number elements will be cast to NaN
    
    if class(DHTable) == "string"
        DHTable = double(DHTable);
    end

    if ~isa(DHTable,'sym') && any(isnan(DHTable), 'all')
        disp("Wrong argument type!")
        disp("Only numeric or sym matrix are supported.")
        return
    end

    % Handle more than one argument
    switch(nargin)
        case 1
            % Build the whole Forward Kinematics
            index = n;
        case 2
            % Up to the specified joint
            index = varargin{2};
        case 3
            % Whole FK, with pre- and post- transformations
            Tb0 = varargin{2};
            Tne = varargin{3};
        case 4
            % With pre- and post- transformations, up to the specified
            % joint
            Tb0 = varargin{2};
            Tne = varargin{3};
            index = varargin{4};
        otherwise
           disp("Too many arguments provided!")
           return
    end
    
    % Check subsequent arguments
    if isnan(index) % TODO: check 'overflow'
        disp("Wrong argument type!")
        disp("Specify number of joint with a number.")
        return
    end

    % Initial homogeneous matrix
    H = eye(4);

    % Handle pre-moltiplication matrix
    if exist("Tb0",'var')
        if any(isnan(Tb0,'all'))
            disp("Tb0 matrix must be numeric!")
            return
        end
        H = Tb0;
    end
    
    % Build Jacobian joint by joint
    for i = [1:index]
       H = H * FKRow(DHTable(i,:));
    end

    % Handle post-moltiplication matrix
    if exist("Tne",'var')
        if any(isnan(Tne,'all'))
            disp("Tne matrix must be numeric!")
            return
        end
        H = H * Tne;
    end
       
    FKmat = H;
end