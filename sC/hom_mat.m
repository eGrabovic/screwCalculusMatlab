function out = hom_mat(rotm, disp)
% HOM_MAT(rotm, disp) returns the homogeneous rototranslational matrix OUT
%   related to the rotation matrix ROTM and the displacement DISP.
%
%   Allows arguments as ADvar instances.

    % Manage ADvars
    if isa(rotm, 'ADvar') || isa(disp, 'ADvar')
        b = ADvar([0 0 0 1], [0 0 0 0]);
        out = [rotm disp; b];

    else % standard case
        % Check to make sure the dimensions of the arguments make sense
        assert(isvector(pos), "pos is not a vector");
        assert(ismatrix(Rot), "Rot is not a matrix");

        n = max(size(pos));
        
        assert(n == 3, ...
            "Array size is not acceptable; it must consist of 3 elements");
        % TODO: substitute with isMatRot when available
        assert(all(sizes(Rot) == [n, n]), ...
            "Matrix and vector sizes are not consistent");

        out = [[rotm,       disp];
               [zeros(1,3),    1]];
    end
            
end