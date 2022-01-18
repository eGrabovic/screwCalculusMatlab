function eqs = DynamicEquations(DHtable, CGtable, Masslist, Tensortable, g0, q, qp, v, vp)
		
	B = Inertia(DHtable, CGtable, Masslist, Tensortable);
	C = InertiaToCoriolis(B, q, qp);
	G = Gravitational(DHtable, CGtable, Masslist, g0);
	
	eqs = B * vp + C * v + G;
end