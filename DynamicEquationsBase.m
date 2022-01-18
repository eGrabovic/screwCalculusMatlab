function eqs = DynamicEquationsBase(DHtable, CGtable, Masslist, Tensortable, gb, q, qp, v, vp, Tb0, Tne)
		
	B = InertiaBase(DHtable, CGtable, Masslist, Tensortable, {Tb0, Tne});
	C = InertiaToCoriolis(B, q, qp);
	G = GravitationalBase(DHtable, CGtable, Masslist, gb, {Tb0, Tne});
		
	eqs = B.vp + C.v + G;
end