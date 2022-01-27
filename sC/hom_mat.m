function out = hom_mat(rotm,disp)
% HOM_MAT(rotm, disp) returns the homogeneous rototranslational matrix OUT
%   related to the rotation matrix ROTM and the displacement DISP.

if isa(rotm, 'ADvar') || isa(disp, 'ADvar')
    b = ADvar([0 0 0 1], [0 0 0 0]);
else
b = [0 0 0 1];
end

out = [rotm disp; b];
            
end