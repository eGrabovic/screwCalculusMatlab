function Gst = hom_mat(screwObj, rotm, disp)
% HOM_MAT(rotm, disp) returns the homogeneous rototranslational matrix OUT
%   related to the rotation matrix ROTM and the displacement DISP.

    if all(size(rotm) ~= [3,3]) || all(size(disp) ~= [3,1])
        error(['Rotation matrix must have 3x3 dimension and the ' ...
                'displacement vector must be a 3x1 array'])
    end

    Gst = [[rotm,           disp];
           [zeros(1, 3),    1]];
    
end