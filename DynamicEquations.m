function eqs = DynamicEquations(DHtable, CGtable, Masslist, Tensortable, ...
                                g0, q, qp, v, vp)
% DYNAMICEQUATIONS(DHtable, CGtable, Masslist, Tensortable, g0, q, qp, v, vp)
%   computes the LHS B(q)vp + C(q,qp)v + G(q) [TODO: how to insert formula?]
%   for the supplied:
%   - DHTABLE (Denavit-Hartemberg table)
%   - CGTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - TENSORTABLE (table with link inertia tensors)
%   - G0 (componenents of gravity in {0})
%   - Q (configuration array)
%   - QP (joint velocity)
%   - V (Slotine-Li velocity)
%   - VP (Slotine-Li acceleration)
%
% For the classical version simply set: v ->qp , vp -> qpp

		
	B = Inertia(DHtable, CGtable, Masslist, Tensortable);
	C = InertiaToCoriolis(B, q, qp);
	G = Gravitational(DHtable, CGtable, Masslist, g0);
	
	eqs = B * vp + C * v + G;
end