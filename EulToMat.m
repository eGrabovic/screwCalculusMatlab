function mat = EulToMat(EU_angles, type, deg_bool)
% EULTOMAT(EU_angles, type, deg_bool) returns the Rotation matrix MAT 
%   related to the specified angles EU_ANGLES and parametrization TYPE.
%   Valid parametrizations consist of a combination of 'X', 'Y', 'Z'
%   (uppercase).
%   If a boolean DEG_BOOL is specified as 'true', input angles will be
%   considered as in degrees.
%   Note that the output matrix is always with angles in radiants.
%
%   Examples:
%       EulToMat([0,0,pi/4], 'XYZ') or EulToMat([0,0,pi/4], "XYZ")
%   -->  [0.7071   -0.7071   0
%        0.7071    0.7071   0
%        0         0        1.0000]
%
%   Input
%       EU_angles:  1x3 numerical array of angles
%       type:       parametrization 'code', chars or string
%       deg_bool:   boolean to indicate that the angles are in degrees
%   Output
%       mat:        3x3 Rotation matrix
%

    assert(all(size(EU_angles) == [1 3]), ...
           "Invalid array of angles specified");

    assert(class(type) == "string" || class(type) == "char", ...
           "'type' argument is not a string or a combination of chars");

    if(class(type) == "char") % handles character arrays
        input_comb = string(type);
    else
        input_comb = type;
    end

    % Build all possible parametrization "codes"
    set = 'XYZ';
    valid_combs = allcomb(set, set, set);

    % Check if a valid Euler 'code' has been specified
    assert(ismember(input_comb, valid_combs), ...
           "Specified combination is not valid")

    if(nargin == 2) % nothing specified
        is_deg = false;
    else
        is_deg = deg_bool;
    end

    % Convert it back to char array for ease of use
    input_comb = char(input_comb);
    mat = eye(3);

    % Build the actual matrix
    for i=[1:size(input_comb,2)]
        switch(input_comb(i))
            case {'Z', 'z'}
                mat = mat * RotZ(EU_angles(i), is_deg);
            case {'X', 'x'}
                mat = mat * RotX(EU_angles(i), is_deg);
            case {'Y', 'y'}
                mat = mat * RotY(EU_angles(i), is_deg);
            otherwise
                disp("Error while parsing the specified combination");
        end
    end
end