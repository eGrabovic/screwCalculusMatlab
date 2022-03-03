function or_err = FramesToVect(R_act, R_des)
% FRAMESTOVECT(R_act, R_des) returns the orientation error OR_ERR of frame
%   R_DES (desired) w.r.t. R_ACT (actual) in the spatial frame where both 
%   R_ACT and R_DES are expressed.
%
%   Examples: TODO
%
%   Input
%       R_act:  3x3 rotation matrix
%       R_des:  3x3 rotation matrix
%   Ouutput
%       or_err: 1x3 orientation error

    % TODO: how to check if it's a rotation matrix?
    assert(all(size(R_act) == [3,3]) || all(size(R_curr) == [3,3]),...
           "Supplied matrices must be 3x3 Rotation matrices");
    
    h_normal = Hat(R_act(:, 1));    
    h_slide = Hat(R_act(:, 2));   
    h_approach = Hat(R_act(:, 3));
    
    normal_des = R_des(:, 1); 
    slide_des = R_des(:, 2); 
    approach_des = R_des(:, 3);
    
    or_err = (1/2) * (h_normal * normal_des + ...
                      h_slide * slide_des + ... 
                      h_approach * approach_des);
end