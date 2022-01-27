function bjac = QuatToBodyJac(quat)
% QUATTOBODYJAC(quat) returns the Body Jacobian BJAC corresponding to the 
%   unit quaternion QUAT (1-vector)
%   If q = (q0, q1, q2, q3) --> w_b = bjac * dot{q}.
% TODO: what is "1-vector"? 1 element?also, hoiw to add formulas?

    b0 = quat(1);
    b1 = quat(2);
    b2 = quat(3);
    b3 = quat(4);
    
    bjac = 2 * [[ -b1,  b0,   b3,  -b2 ]
		        [ -b2, -b3,   b0,   b1 ]
		        [ -b3,  b2,  -b1,   b0 ]];
end
	