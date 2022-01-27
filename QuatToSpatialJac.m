function sjac = QuatToSpatialJac(quat)
% QUATTOSPATIALJAC(quat) returns the Spatial Jacobian SJAC corresponding to
%   the unit quaternion QUAT (1-vector).
%   If q = (q0, q1, q2, q3), w_s = Js dot{q}."
%   TODO: 1-vector? formulas?

    b0 = quat(1);
    b1 = quat(2);
    b2 = quat(3);
    b3 = quat(4);
    
    sjac = 2 * [[ -b1,  b0, -b3,  b2 ]
	    	    [ -b2,  b3,  b0, -b1 ]
		        [ -b3, -b2,  b1,  b0 ]];
end