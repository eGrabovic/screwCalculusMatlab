%% Graphical Primitives
% 
% CreateGround({xmin_, xmax_, dx_}, {ymin_, ymax_, dy_}, colorplane_, colorframe_) 
% 	%Module({x, y, planepts, plane, e1s, e2s, e3s, Orig, frameS, scene},
% 		
% 		planepts = Table({x, y, 0}, {x, xmin, xmax, dx}, {y, ymin, ymax, dy});
% 
%         plane = Map( Graphics3D, {Map({colorplane, Line(#)} &, planepts), Map({colorplane, Line(#)} &, Transpose(planepts))});
%         
%         e1s = {5, 0, 0}; e2s = {0, 5, 0}; e3s = {0, 0, 5}; Orig = {0, 0, 0};
%         
%         frameS = Graphics3D({Thick, colorframe, Line({Orig, e1s}), Line({Orig, e2s}), Line({Orig, e3s})});
% 
%         scene = Append(plane, frameS);
%         
%         Return(scene);
% 		
% 	);
% 
% CreateFrame(g_, colorframe_) 
% 	%Module({R, p, x, y, z, scale, frame},
% 		scale = 5.0;
% 		
% 		R = RigidOrientation(g);
% 		p = RigidPosition(g);
% 		
% 		px = p + scale Flatten(TakeColumns(R, {1}));
% 		py = p + scale Flatten(TakeColumns(R, {2}));
% 		pz = p + scale Flatten(TakeColumns(R, {3}));
% 		
% 		frame = Graphics3D({Thick, colorframe, Line({p, px}), Line({p, py}), Line({p, pz})});
% 	
% 	    Return(frame);
% 	);
% 	
% CreateJoint(g_, {height_, radius_}) 
% 	%Module({R, p, z, bottom, top, cylinder},
% 		
% 		R = RigidOrientation(g);
% 		p = RigidPosition(g);
% 		z = Flatten(TakeColumns(R, {3}));
% 		
% 		bottom = p - (height/2) z;
% 		top    = p + (height/2) z;
% 		
% 		cylinder = Graphics3D({
% 			                   Green, Opacity(0.25), Specularity(White, 30),
% 		                       Cylinder({bottom, top}, radius)
% 		                       }, Boxed -> False);
% 		
% 		Return(cylinder);	
% 	);
% 	
% CreateLink(g_, {length_, radius_}) 
% 	%Module({R, p, x, start, end, cylinder},
% 		
% 		R = RigidOrientation(g);
% 		p = RigidPosition(g);
% 		x = Flatten(TakeColumns(R, {1}));
% 		
% 		start  = p - length x;
% 		end    = p;
% 		
% 		cylinder = Graphics3D({
% 			                   Orange, Opacity(0.25), Specularity(White, 30),
% 		                       Cylinder({start, end}, radius)
% 		                       }, Boxed -> False);
% 		
% 		Return(cylinder);	
% 	);
% 	
% CreateRobot(DHtable_, {gpb_, gb0_}, q0_, a0_:20) 
% 	%Module({Tpb, Tb0, i, n, rlink, rjoint, hjoint, robot, res},
% 		Tpb = gpb;
% 		Tb0 = gb0;
% 		n   = Length(DHtable@@q0);
% 		rlink  = 2.5;
% 		rjoint = 5.0; 
% 		hjoint = 10.0;
% 		     % a0 = 20; *)
% 		
% 		% ruleq = Thread(Rule(q, q0)); *)
% 		
% 		link0  = {CreateLink(Tpb.Tb0.HomogeneousRotY(-Pi/2), {a0, rlink})};
% 		joint1 = {CreateJoint(Tpb, {hjoint, rjoint})};
% 		
% 		linklist = Table( 
% 			             CreateLink( 
% 			             	   DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)}, i 
% 			             	          ), 
% 			             		        { (DHtable@@q0)((i, 1)), rlink } 
% 			                       ),
% 			             {i, 1, n}
% 			             );
% 			             
% 		jointlist = Table( 
% 			             CreateJoint( 
% 			             	   DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)}, i-1 
% 			             	          ), 
% 			             		        { hjoint, rjoint } 
% 			                       ),
% 			             {i, 2, n}
% 			             );	             
% 		frameEE = {CreateFrame( DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)} ), Magenta )};
% 			             
% 		res = Join(link0, joint1, linklist, jointlist, frameEE);	             
% 		
% 		Return(res);
% 	);
% 
% BBox(a_, b_, c_) 
% 	%Module({list},
% 		p1 = {a, b, -c, 1}; p2 = {-a, b, -c, 1}; p3 = {-a, -b, -c, 1}; p4 = {a, -b, -c, 1};
% 		p5 = {a, b,  c, 1}; p6 = {-a, b,  c, 1}; p7 = {-a, -b,  c, 1}; p8 = {a, -b,  c, 1};
% 	    
% 	    list1 = {p1, p2, p3, p4, p1}//Transpose;
% 	    list2 = {p5, p6, p7, p8, p5}//Transpose;
% 	    list3 = {p1, p5}//Transpose;
% 	    list4 = {p2, p6}//Transpose;
% 	    list5 = {p3, p7}//Transpose;
% 	    list6 = {p4, p8}//Transpose;
% 	    
% 	    list = {list1, list2, list3, list4, list5, list6};
% 	    
% 	    Return(list);
% 	
% 	);
% 
% EEllipsoid(a_, b_, c_) 
% 	%Module({nu, nv},
% 		
% 		nu = 20.; nv = 20.;
% 		
% 		x(u_, v_)  a Cos(v) Cos(u);
% 		y(u_, v_)  b Cos(v) Sin(u);
% 		z(u_, v_)  c Sin(v);
% 		
% 		r(u_, v_)  {x(u,v), y(u,v), z(u,v), 1};
% 		
% 		list1 = Table(r(u,v), {u, 0, 2. Pi, 2. Pi/nu}, {v, -Pi/2., Pi/2., Pi/nv});
% 		list2 = Transpose(list1);
% 		list1 = Transpose /@ list1;
% 		list2 = Transpose /@ list2;
% 		Return( Join(list1, list2) );
% 		
% 	);
% 	
% CreateObject(g_, objectlist_, color_) 
% 	%Module({points, lines, res},
% 		
% 		points = Transpose(Drop(Dot(g,#),-1)) & /@ objectlist;
% 		
% 		lines = Line /@ points;
% 		
% 		res = Graphics3D({Thin, color, lines});
% 		
% 		Return(res);
% 		
% 	);
% 	
% FindRedundantAnglesFromOrigin(origin_, DHtable_, {Tb0_, Tne_}, q_, qguess_) 
% 	%Module({g, p, f, L, gcons1, gcons2, gcons, cons, initguess, res},
% 		
% 		g = DHFKine(DHtable, {Tb0, Tne});
% 		p = RigidPosition(g);
% 		
% 		f      = {(1/2) q.q}; 
% 		% f      = {Max(q)};    *)
% 		
% 		gcons1  = origin - p;
% 		gcons2  = {};
% 		gcons   = Join(gcons1, gcons2);
% 		cons     = Thread(Equal(gcons, 0));
% 		
% 		L = Join(f, cons);
% 		
% 		initguess = Thread(List(q, qguess));
% 		
% 		res = q /. FindMinimum(L, initguess )((2)); 
% 		% res = q /. NMinimize(L, q )((2)); *)
% 
%         Return(res);		
% 	);

