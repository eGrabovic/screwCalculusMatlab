function eqs = DynamicEquationsBase(DHtable, CGtable, Masslist, Tensortable, ...
                                    gb, q, qp, v, vp, Tb0, Tne)
% DYNAMICEQUATIONSBASE(DHtable, CGtable, Masslist, Tensortable, gb, q, qp, v, vp, Tb0, Tne)
%   computes the LHS B(q)vp + C(q,qp)v + G(q) [TODO: how to insert formula?]
%   for the supplied:
%   - DHTABLE (Denavit-Hartemberg table)
%   - CGTABLE (pici table)
%   - MASSLIST (list of link masses)
%   - TENSORTABLE (table with link inertia tensors)
%   - GB (componenents of gravity in {B})
%   - Q (configuration array)
%   - QP (joint velocity)
%   - V (Slotine-Li velocity)
%   - VP (Slotine-Li acceleration)
%   - [TB0, TNE] initial [B to S0] and final [Sn to E] offset transformations
%
%   For the classical version simply set: v ->qp , vp -> qpp.

		
	B = InertiaBase(DHtable, CGtable, Masslist, Tensortable, {Tb0, Tne});
	C = InertiaToCoriolis(B, q, qp);
	G = GravitationalBase(DHtable, CGtable, Masslist, gb, {Tb0, Tne});
		
	eqs = B.vp + C.v + G;
end