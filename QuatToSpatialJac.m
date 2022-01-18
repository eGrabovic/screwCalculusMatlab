function sjac = QuatToSpatialJac(quat)
    b0 = quat(1);
    b1 = quat(2);
    b2 = quat(3);
    b3 = quat(4);
    
    sjac = 2 * [[ -b1,  b0, -b3,  b2 ]
	    	    [ -b2,  b3,  b0, -b1 ]
		        [ -b3, -b2,  b1,  b0 ]];
end