%%
% Rodriguez parameters gamma = r tan(theta/2) *)
RodriguezToMat(gamma1_, gamma2_, gamma3_)  
	%Module( {gamma, hatgamma, modulusgammasquared},
	
			gamma = {gamma1, gamma2, gamma3};
			hatgamma = Hat(gamma);
			modulusgammasquared = gamma.gamma;
			
			IdentityMatrix(3) + 2/(1 + modulusgammasquared) (hatgamma + hatgamma.hatgamma)
	      );
	      
RodriguezToSpatialJac(gamma1_, gamma2_, gamma3_)  
	%Module( {gamma, modulusgammasquared},
		gamma = {gamma1, gamma2, gamma3};
		modulusgammasquared = gamma.gamma;
		2/(1 + modulusgammasquared) { {      1,   -gamma3,    gamma2},
									  { gamma3,         1,   -gamma1},
									  {-gamma2,    gamma1,         1}
									}  	
	
	);
	
RodriguezToBodyJac(gamma1_, gamma2_, gamma3_)  
	%Module( {gamma, modulusgammasquared},
		gamma = {gamma1, gamma2, gamma3};
		modulusgammasquared = gamma.gamma;
		2/(1 + modulusgammasquared) { {      1,   gamma3,    -gamma2},
									  { -gamma3,         1,   gamma1},
									  {  gamma2,    -gamma1,         1}
									}  	
	
	);	

%
%  *  Spatial and Body Jacobians relating Spatial or Body fixed angular velocity components to
%  *  angle derivatives in various parametrizations
%  *)

EulZXZToSpatialJac(phi_, theta_, psi_) 
	%Module({},
		{{ 0, Cos(phi),  Sin(phi) Sin(theta) },
		 { 0, Sin(phi), -Cos(phi) Sin(theta) },
		 { 1,        0,           Cos(theta) }
		}
		
	);
	
EulZXZToBodyJac(phi_, theta_, psi_) 
	%Module({},
		{{ Sin(theta) Sin(psi),  Cos(psi), 0 },
		 { Sin(theta) Cos(psi), -Sin(psi), 0 },
		 {          Cos(theta),         0, 1 }
		}
		
	);	
	
EulZYZToSpatialJac(phi_, theta_, psi_) 
	%Module({},
		{{ 0, -Sin(phi), Cos(phi) Sin(theta)  },
		 { 0,  Cos(phi), Sin(phi) Sin(theta)  },
		 { 1,        0,           Cos(theta)  }
		}
		
	);	

EulZYZToBodyJac(phi_, theta_, psi_) 
	%Module({},
		{{ -Cos(psi) Sin(theta), Sin(psi), 0 },
		 {  Sin(psi) Sin(theta), Cos(psi), 0 }
		 {           Cos(theta),        0, 1 }
		}
		
	);	
	
EulZYXToSpatialJac(phi_, theta_, psi_) 
	%Module({},
		{{ 0, -Sin(phi), Cos(theta) Cos(phi)  },
		 { 0,  Cos(phi), Cos(theta) Sin(phi)  },
		 { 1,        0,          -Sin(theta)  }
		}
		
	);		
	
EulZYXToBodyJac(phi_, theta_, psi_) 
	%Module({},
		{{          -Sin(theta),         0,  1 }         
		 {  Cos(theta) Sin(psi),  Cos(psi),  0 }
		 {  Cos(theta) Cos(psi), -Sin(psi),  0 }
		}
		
	);		
	
QuatToSpatialJac( { b0_, b1_, b2_, b3_} )  
	%Module({},
	
	2 * { { -b1,  b0, -b3,  b2  },
		  { -b2,  b3,  b0, -b1 },
		  { -b3, -b2,  b1,  b0 }
	    }
	
	);	*)

QuatToSpatialJac( b0_, b1_, b2_, b3_ )  
	%Module({},
	
	2 * { { -b1,  b0, -b3,  b2  },
		  { -b2,  b3,  b0, -b1 },
		  { -b3, -b2,  b1,  b0 }
	    }
	
	);
	
% QuatToBodyJac( { b0_, b1_, b2_, b3_} )  
	%Module({},
	
	2 * { { -b1,  b0,  b3, -b2  },
		  { -b2, -b3,  b0,  b1 },
		  { -b3,  b2, -b1,  b0 }
	    }
	
	);		*)

QuatToBodyJac( b0_, b1_, b2_, b3_ )  
	%Module({},
	
	2 * { { -b1,  b0,  b3, -b2  },
		  { -b2, -b3,  b0,  b1 },
		  { -b3,  b2, -b1,  b0 }
	    }
	);		
%%









WrenchAxis(xi_?VectorQ)   Null;

% Find the magnitude of a twist *)
TwistMagnitude(xi_?VectorQ)   
  %Module(
    {v, w},
    {v, w} = Partition(xi, 3);
    if ((MatchQ(w,{0,0,0}) || MatchQ(w, {{0},{0},{0}})), 
      Sqrt(v.v), Sqrt(w.w))
  );
WrenchMagnitude(xi_?VectorQ)   Null;
 

%
 * Adjoint calculation
 *
 * The adjoint matrix maps twist vectors to twist vectors.
 *
 *)

% Adjoint matrix calculation *)
RigidAdjoint(g_?MatrixQ)   
  %Module(
    {R = RigidOrientation(g), p = RigidPosition(g)},
    
    BlockMatrix({{R,       AxisToSkew(p) . R},
                 {ZeroMatrix(3,3),  R}
                 }
    )
    
  );
  
% Inverse adjoint matrix calculation *)
InverseRigidAdjoint(g_?MatrixQ)   
  %Module(
    {RT = Transpose(RigidOrientation(g)), p = RigidPosition(g)},
    
    BlockMatrix({{RT,       -RT.AxisToSkew(p)},
                 {ZeroMatrix(3,3),  RT}
                 }
    )
    
  );  
  
 % Calculation of the error twist xi_err and error angle theta_err such that:
    gStart.TwistExp(xi_err, theta_err) = gEnd *)
    
 LocalTranTwist(gStart_?MatrixQ, gEnd_?MatrixQ) 
 %Module(
     {gError, xi, theta, i},
     
     gError = Inverse(gStart).gEnd;
     
     if (gError == IdentityMatrix(4),
       Return({Table(0,{i,6}),0});
     );
     
     {xi,theta} = RigidTwist(gError);
     
     Return({xi,theta});
);

% Calculation of the error twist xi_err and error angle theta_err such that:
    TwistExp(xi_err, theta_err).gStart = gEnd *)
  
 GlobalTranTwist(gStart_?MatrixQ, gEnd_?MatrixQ) 
 %Module(
     {gError, xi, theta, Adg, i},
     
     gError = Inverse(gStart).gEnd;
     
     if (gError == IdentityMatrix(4),
       Return({Table(0,{i,6}),0});
     );
     
     {xi,theta} = RigidTwist(gError);
     
     Adg = RigidAdjoint(gStart);
      xi = Adg.xi;
       
     Return({xi,theta});
);
  
%
 *  Quaternions to matrix and matrix to quaternions conversion utilities
 *
 *
 *)



% ::Section:: *)
%Quaternion Library:*)
%   q((1))        -->      q0 = Cos(th/2)*)
%   q((2,3,4))) -->  qvec = Sin(th/2) nhat*)


% ::Subsection:: *)
%ROTATIONS and quaternions:*)
%q*v*qbar  with v= (0,vec) a "pure" quaternion*)
%gives EACH COLUMN of a rotation matrix:*)





% ::Subsection:: *)
%Warning: these are the reverse of conventional argument order:*)
%change order once all is consistent to (angle, nhat).*)


MakeQuat(n_List:{0,0,1},angle_:0)  
     %Module({c = Cos(angle/2), s = Sin(angle/2)},
     {c, n((1))*s, n((2))*s, n((3))*s}//N)


MakeRot(n_List:{0,0,1},angle_:0)  
%Module({c = Cos(angle)//N, s = Sin(angle)//N, cm},
cm = 1 - c;
{{c + cm*n((1))^2, cm*n((2))*n((1)) - s*n((3)), cm*n((3))*n((1)) + s*n((2))},
{cm*n((1))*n((2)) + s*n((3)), c + cm*n((2))^2,  cm*n((3))*n((2)) - s*n((1))},
{ cm*n((1))*n((3)) - s*n((2)), 
              cm*n((2))*n((3)) + s*n((1)), c + cm*n((3))^2}}//N)


RotToQuatSym(mat_List)  
   %Module({q0,q1,q2,q3,trace,s,t1,t2,t3},
         trace = Sum(mat((i,i)),{i,1,3});
             s = Sqrt(trace + 1.0);
               q0 = s/2;
               s = 1/(2 s);
               q1 = (mat((3,2)) - mat((2,3)))*s;
               q2 = (mat((1,3)) - mat((3,1)))*s;
               q3 = (mat((2,1)) - mat((1,2)))*s;
                {q0,q1,q2,q3})

%






% ::Subsection:: *)
%Distances in quaternion arcs:*)


QLength(curveOnSphere_List)  
  Sum(ArcCos(curveOnSphere((i)) . curveOnSphere((i+1))),
        {i,1,Length(curveOnSphere)-1})


% ::Subsection:: *)
%Arrays of Quaternion frames need these:  *)


ForceCloseQFrames(qfrms_List)  
  %Module({lastfrm = qfrms((1)), thisfrm},
    Table(thisfrm = qfrms((i));
          lastfrm =
             if (thisfrm . lastfrm >= 0, thisfrm, - thisfrm),
             {i,1,Length(qfrms)}))


 ForceClose2DQFrames(qrows_List)  
  %Module({lastrow = qrows((1)), thisrow},
    Table(thisrow = ForceCloseQFrames(qrows((i)));
          lastrow = 
             if (thisrow((1)) . lastrow((1)) >= 0, thisrow, - thisrow),
              {i,1,Length(qrows)}))
              


ForceClose3DQFrames(qplanes_List) %Module({lastplane=qplanes((1)),thisplane},Table(thisplane=ForceClose2DQFrames(qplanes((i)));
lastplane=if (thisplane((1)).lastplane((1))>=0,thisplane,-thisplane),{i,1,Length(qplanes)}))

*)



% ::Subsection:: *)
%4D rotation double quaternions*)


QQTo4DRot(p_List,q_List)  
  %Module({q0= q((1)), q1=q((2)), q2 = q((3)), q3 = q((4)),
         p0= p((1)), p1=p((2)), p2 = p((3)), p3 = p((4))},
{{p0*q0 + p1*q1 + p2*q2 + p3*q3, p1*q0 - p0*q1 - p3*q2 + p2*q3, 
            p2*q0 + p3*q1 - p0*q2 - p1*q3, p3*q0 - p2*q1 + p1*q2 - p0*q3}, 
 {(-(p1*q0) + p0*q1 - p3*q2 + p2*q3), (p0*q0 + p1*q1 - p2*q2 - p3*q3), 
     (-(p3*q0) + p2*q1 + p1*q2 - p0*q3), (p2*q0 + p3*q1 + p0*q2 + p1*q3)}, 
 {(-(p2*q0) + p3*q1 + p0*q2 - p1*q3), (p3*q0 + p2*q1 + p1*q2 + p0*q3),
     (p0*q0 - p1*q1 + p2*q2 - p3*q3), (-(p1*q0) - p0*q1 + p3*q2 + p2*q3)}, 
 {(-(p3*q0) - p2*q1 + p1*q2 + p0*q3), (-(p2*q0) + p3*q1 - p0*q2 + p1*q3), 
     (p1*q0 + p0*q1 + p3*q2 + p2*q3), (p0*q0 - p1*q1 - p2*q2 + p3*q3)}})
 


% ::Section:: *)
%Spline Library is here*)


Lerp(p0_List,p1_List,t_)   (1-t)*p0 + t*p1


Slerp(p0_List,p1_List,t_)   %Module({costh = p0 . p1//Chop, th, sinth},
  if (costh > 0.0,costh = Chop(costh -1.) +1.);
  if (costh < 0.0,costh = Chop(costh +1.) -1.);
 th = N(ArcCos(costh));
 sinth = N(Sin(th));
if (sinth == 0, 
 (1-t)*p0 + t*p1,
 (Sin(th*(1-t))/sinth)*p0 +(Sin(th*t)/sinth)*p1 ))


% ::Section:: *)
%Spline and Spherical Splines*)


CRspline(p0_,p1_,p2_,p3_,t_)  
	Lerp(Lerp(Lerp(p0,p1,t+1),Lerp(p1,p2,t),(t+1)/2.),
              Lerp(Lerp(p1,p2,t),Lerp(p2,p3,t-1),           t/2.),         t)


BZspline(p0_,p1_,p2_,p3_,t_)  
	Lerp(Lerp(Lerp(p0,p1,t),Lerp(p1,p2,t),t),
              Lerp(Lerp(p1,p2,t),Lerp(p2,p3,t),t),t)


UBspline(p0_,p1_,p2_,p3_,t_)  
	Lerp(Lerp(Lerp(p0,p1,(t+2)/3.),Lerp(p1,p2,(t+1)/3.),(t+1)/2.),
              Lerp(Lerp(p1,p2,(t+1)/3.),Lerp(p2,p3,t/3.),                       t/2.),t)


SCRspline(p0_,p1_,p2_,p3_,t_)  
	Slerp(Slerp(Slerp(p0,p1,t+1),Slerp(p1,p2,t),(t+1)/2.),
              Slerp(Slerp(p1,p2,t),Slerp(p2,p3,t-1),           t/2.),         t)


SBZspline(p0_,p1_,p2_,p3_,t_)  
	Slerp(Slerp(Slerp(p0,p1,t),Slerp(p1,p2,t),t),
              Slerp(Slerp(p1,p2,t),Slerp(p2,p3,t),t),t)


SUBspline(p0_,p1_,p2_,p3_,t_)  
	Slerp(Slerp(Slerp(p0,p1,(t+2)/3.),Slerp(p1,p2,(t+1)/3.),(t+1)/2.),
              Slerp(Slerp(p1,p2,(t+1)/3.),Slerp(p2,p3,t/3.),                       t/2.),t)


squad(x0_,x1_,x2_,x3_,t_)  (1-2t (1 - t))*((1-t) x0 + t x3) + 2t*(1-t)((1-t)x1 + t x2)

%
 *  Calculations of the forward kinematic map and
 *  the Spatial and Body Jacobians
 *
 *)
 
 % Gives Xi 6 vector given a point on axis and axis unit vector for a Revolute Joint *)
RevoluteTwist(q_,w_)  Flatten({Cross(q,w),w});

% Gives Xi 6 vector given a point on axis and axis unit vector for a Prismatic Joint *)
PrismaticTwist(q_,w_)  Flatten({w, {0,0,0}});

% Gives the homogeneous matrix *)
ForwardKinematics(args__, gst0_)  
  %Module({ g, i,
      argList = {args},		% turn arguments into a real list *)
      n = Length({args})	% decide on the number of joints *)
    },

    % Initialize the transformation matrix *)
    g = TwistExp(argList((1,1)), argList((1,2)));   

    % Build up the Jacobian joint by joint *)
    For(i = 2, i <= n, i++,
      % Update the transformation matrix *)
      g = g . TwistExp(argList((i,1)), argList((i,2)));
    );      

    % Finish by multiplying by the initial tool configuration *)
    g . gst0
  );
			

% Construct the Spatial Jacobian for a robot with any no. of links *)
SpatialJacobian(args__, gst0_)   
  %Module(
    { i, xi, Js, g,
      argList = {args},		% turn arguments into a real list *)
      n = Length({args})	% decide on the number of joints *)
    },

    % First initialize the Jacobian and compute the first column *)
    Js = {argList((1,1))};
    g = TwistExp(argList((1,1)), argList((1,2)));   

    % Build up the Jacobian joint by joint *)
    For(i = 2, i <= n, i++,
      % Compute this column of the Jacobian and append it to Js *)
      xi = RigidAdjoint(g) . argList((i,1));
      Js = Join(Js, {xi});      

      % Update the transformation matrix *)
      g = g . TwistExp(argList((i,1)), argList((i,2)));
    );      

    % Return the Jacobian *)
    Transpose(Js)
  );
			
% Construct the Body Jacobian for a robot with any no. of links *)			
BodyJacobian(args__, gst0_)   
  %Module(
    { i, xi, Jb, g,
      argList = {args},		% turn arguments into a real list *)
      n = Length({args})	% decide on the number of joints *)
    },

    % Initialize the Jacobian and the transformation matrix *)
    Jb = {};
    g = gst0;

    % Build up the Jacobian joint by joint *)
    For(i = n, i >= 1, i--,
      % Compute this column of the Jacobian and prepend it to Jb *)
      xi = RigidAdjoint(RigidInverse(g)) . argList((i,1));
      Jb = Join({xi}, Jb);      

      % Update the transformation matrix *)
      g = TwistExp(argList((i,1)), argList((i,2))) . g;
    );      

    % Return the Jacobian *)
    Transpose(Jb)
  );
  
% Definitions for Denavit-Hartenberg convention *)

HomogeneousRotX(alpha_)  
	%Module({},
		R = RotX(alpha);
		d = {0,0,0};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);
  
HomogeneousRotY(alpha_)  
	%Module({},
		R = RotY(alpha);
		d = {0,0,0};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);

HomogeneousRotZ(alpha_)  
	%Module({},
		R = RotZ(alpha);
		d = {0,0,0};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);
	
HomogeneousTranslX(alpha_)  
	%Module({},
		R = IdentityMatrix(3);
		d = {alpha,0,0};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);	
	
HomogeneousTranslY(alpha_)  
	%Module({},
		R = IdentityMatrix(3);
		d = {0, alpha, 0};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);	

HomogeneousTranslZ(alpha_)  
	%Module({},
		R = IdentityMatrix(3);
		d = {0, 0, alpha};
		H = RPToHomogeneous(R, d);
		Return(H);  
	);
	
DH(pars_List)  
	%Module({   a = pars((1)), 
		   alpha = pars((2)),
		       d = pars((3)),
		   theta = pars((4)),
	   jointtype = pars((5)),
		   TZ, RZ, TX, RX},
		       TZ = HomogeneousTranslZ(d); 
		       RZ = HomogeneousRotZ(theta);
		       TX = HomogeneousTranslX(a);
		       RX = HomogeneousRotX(alpha);
		        H = TZ.RZ.TX.RX;
		        Return(H);
	);

%{ 
DHFKine(DHtable_)  
  %Module(
    { H, i, n = Length(DHtable)	% decide on the number of joints *)
    },

    % Initialize homogeneous matrix *)
    H = IdentityMatrix(4);

    % Build up the Jacobian joint by joint *)
    For(i = 1, i <= n, i++,
      % Update the transformation matrix *)
      H = H . DH( DHtable((i)) );
    );      

    % Once finished... *)
    Return(H);
  );
  
DHFKine(DHtable_, j_)  
  %Module(
    { H, i, n = Length(DHtable)	% decide on the number of joints *)
    },

    % Initialize homogeneous matrix *)
    H = IdentityMatrix(4);

    % Build up the Jacobian joint by joint *)
    For(i = 1, i <= j, i++,
      % Update the transformation matrix *)
      H = H . DH( DHtable((i)) );
    );      

    % Once finished... *)
    Return(H);
  );  
  
DHFKine(DHtable_, {Tb0_, Tne_}) 
	%Module(
		{H0n, H},
        H0n = DHFKine(DHtable);
        H   = Tb0.H0n.Tne;
        
        Return(H);		
	);
  
DHFKine(DHtable_, {Tb0_, Tne_}, j_) 
	%Module(
    { H, i, n = Length(DHtable)	% decide on the number of joints *)
    },

    % Initialize homogeneous matrix *)
    H = Tb0;

    % Build up the Jacobian joint by joint *)
    For(i = 1, i <= j, i++,
      % Update the transformation matrix *)
      H = H . DH( DHtable((i)) );
    );      

    % Once finished... *)
    Return(H);
  ); 
%}

DHJacob0(DHtable_) 
	%Module({ i, n = Length(DHtable), He,
		     type, z, r,
		     jv, jo, j, J },
		     
		J  = {};     
		He = DHFKine(DHtable);
		For(i = 1, i <= n, i++,
			
			type = DHtable((i, 5));
			   z = (RigidOrientation(DHFKine(DHtable, i-1)))((All, 3));
			   r = RigidPosition(He) - (RigidPosition(DHFKine(DHtable, i-1)));
			
			if ( type == "P",
				
				% Prismatic Joint *)
				jv = z;
				jo = {0, 0, 0};
				j  = Join(jv, jo) ; ,
			
			    % Revolute Joint *)
			    jv = Hat(z).r;
			    jo = z;
			    j  = Join(jv, jo); 
		      );
		J = Append(J, j);
		
	);
	Return(Transpose(J));
	);

DHJacobBase(DHtable_, Tb0_) 	
	%Module({g, R, zero, Adg, J, Jb},
		R    = RigidOrientation(Tb0);
		zero = {0, 0, 0};
		g    = RPToHomogeneous(R, zero);
		Adg  = RigidAdjoint(g);
		J    = DHJacob0(DHtable);
		Jb   = Adg.J;
		
		Return(Jb);
	);
	
DHJacob0Dyn(DHtable_, k_) 	
	%Module({ i, n = Length(DHtable), Hk,
		     type, z, r, 
             jv, jo, j, Jk, J },
		     
		Jk  = {};     
		Hk = DHFKine(DHtable, k);
		For(i = 1, i <= k, i++,
			
			type = DHtable((i, 5));
			   z = (RigidOrientation(DHFKine(DHtable, i-1)))((All, 3));
			   r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, i-1)));
			
			if ( type == "P",
				
				% Prismatic Joint *)
				jv = z;
				jo = {0, 0, 0};
				j  = Join(jv, jo) ; ,
			
			    % Revolute Joint *)
			    jv = Hat(z).r;
			    jo = z;
			    j  = Join(jv, jo); 
		      );
		Jk = Append(Jk, j);
		
		);
		
	    J = StackCols( Transpose(Jk), ZeroMatrix(6, n-k) );
	    Return(J);
	);
	
DHJacobBaseDyn(DHtable_, {Tb0_, Tne_}, k_) 	
	%Module({ i, n = Length(DHtable), Hk,
		     type, z, r, 
             jv, jo, j, Jk, J },
		     
		Jk  = {};    
		Hk  = DHFKine(DHtable, {Tb0, Tne}, k);
	
		For(i = 1, i <= k, i++,
			
			type = DHtable((i, 5));
			   z = (RigidOrientation(DHFKine(DHtable, {Tb0, Tne}, i-1)))((All, 3));
			   r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, {Tb0, Tne}, i-1)));
			
			if ( type == "P",
				
				% Prismatic Joint *)
				jv = z;
				jo = {0, 0, 0};
				j  = Join(jv, jo) ; ,
			
			    % Revolute Joint *)
			    jv = Hat(z).r;
			    jo = z;
			    j  = Join(jv, jo); 
		      );
		Jk = Append(Jk, j);
		
		);
		
	    J = StackCols( Transpose(Jk), ZeroMatrix(6, n-k) );
	    Return(J);
	);	

CGJacob0Dyn(DHtable_, CGtable_, k_) 
	%Module({R0k, pkck, M, DHJacob, CGJacob},
		R0k  = RigidOrientation( DHFKine(DHtable, k) );
		pkck = CGtable((k));
		M = BlockMatrix( {
			               {   IdentityMatrix(3),     -Hat(R0k.pkck) },
						   {     ZeroMatrix(3,3),  IdentityMatrix(3) }
						  }
			           );
		DHJacob = DHJacob0Dyn(DHtable, k);
		CGJacob = M.DHJacob;
		Return(CGJacob);	           
	);

CGJacobBaseDyn(DHtable_, CGtable_, {Tb0_, Tne_}, k_) 
	%Module({R0k, pkck, M, DHJacob, CGJacob},
		R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
		pkck = CGtable((k));
		M = BlockMatrix( {
			               {   IdentityMatrix(3),     -Hat(R0k.pkck) },
						   {     ZeroMatrix(3,3),  IdentityMatrix(3) }
						  }
			           );
		DHJacob = DHJacobBaseDyn(DHtable, {Tb0, Tne}, k);
		CGJacob = M.DHJacob;
		Return(CGJacob);	           
	);

% Calcolo dei termini che formano il regressore *)
TensorDerivative(J_, q_) 	
	%Module({DJ, res },
		DJ  = D(J, {q});
		res =(1/2) (DJ + Transpose(DJ, {1, 3, 2}));
		Return(res);
	);
 
Regressor(DHtable_, q_, qp_, v_, vp_, t_, g0_, k_) 
	%Module({X0p, X1a, X1b, X1ap, X1bp, X1p, X2p, par,
		    W0, W1, W2,
		    Z0, Z1,
		    R0k, R0kT,
		    Jk, Jvk, Jok, JvkT, JokT,
		    Jvkp, JvkpT, Jokp, JokpT,
		    JvTJv,
		    T1, T2, T3, TT, 
		    E1, E2, E3, E4, E5, E6, EE,
		    Y0, Y1, Y2, Y},
		    
		    R0k  = RigidOrientation( DHFKine(DHtable, k) );
		    R0kT = Transpose(R0k);
		    Jk   =  DHJacob0Dyn(DHtable, k);
		    Jvk  =  Jk((1;;3, All));
		    Jok  =  Jk((4;;6, All));
		    
		    JvkT = Transpose(Jvk);
		    JokT = Transpose(Jok);
		    JvTJv = JvkT.Jvk;
		    
		    Jvkp  = TensorDerivative(Jvk, q).qp;
		    JvkpT = Transpose(Jvkp);
		    
		    Jokp  = TensorDerivative(Jok, q).qp;
		    JokpT = Transpose(Jokp); 
		    
		    % termini da d   dT1                                        *)
		    %            -   -                                          *)
		    %            dt  dqp         con T1 = v.B(q).qp             *)
		    
		    
 			X0p = Transpose({(TensorDerivative(JvTJv, q).qp).v  + JvTJv.vp});
 			
 			
 			% X1p = (   JvkpT.Hat(Jok.v) + JvkT.(Hat(Jokp.v) + Hat(Jok.vp))
 				   - JokpT.Hat(Jvk.v)  - JokT.(Hat(Jvkp.v) + Hat(Jvk.vp)) ).R0k +
 				      (JvkT.Hat(Jok.v) - JokT.Hat(Jvk.v)).(D(R0k, t)); *)
 			
 			T1 = SparseArray({{2, 3} -> -1,  {3, 2} ->  1, {3,3} -> 0});
 			T2 = SparseArray({{1, 3} ->  1,  {3, 1} -> -1, {3,3} -> 0});
 			T3 = SparseArray({{1, 2} -> -1,  {2, 1} ->  1, {3,3} -> 0});
 			
 			TT = {T1, T2, T3};
 			TT = Transpose(TT, {2,1,3});
 			
 			X1a  = JokT.R0k.TT.R0kT.Jvk;
 			
 			X1ap = Transpose(Table((TensorDerivative(X1a((All,i,All)), q).qp).v, {i,1,3}));
 			
 			X1b  = -JvkT.R0k.TT.R0kT.Jok;  
 			
 			X1bp = Transpose(Table((TensorDerivative(X1b((All,i,All)), q).qp).v, {i,1,3}));
 			
 			
 			X1p = (X1a + X1b).vp + X1ap + X1bp;
 			
 			E1 = SparseArray({{1, 1} -> 1,  {3, 3} -> 0});
 			E2 = SparseArray({{1, 2} -> -1, {2, 1} -> -1 , {3, 3} -> 0});
 			E3 = SparseArray({{1, 3} -> -1, {3, 1} -> -1 , {3, 3} -> 0});
 			E4 = SparseArray({{2, 2} -> 1,  {3, 3} -> 0});
 			E5 = SparseArray({{2, 3} -> -1, {3, 2} -> -1 , {3, 3} -> 0});
 			E6 = SparseArray({{3, 3} -> 1});
 			EE = {E1, E2, E3, E4, E5, E6};
 			EE = Transpose(EE, {2,1,3});
 			
 			par = JokT.R0k.EE.R0kT.Jok;
 			
 			X2p = Transpose(Table((TensorDerivative(par((All,i,All)), q).qp).v, {i,1,6})) + par.vp;
 		
 			
 			% termini da dT2                                    *)
		    %            -                                      *)
		    %            dq           con T2 = (1/2)v.B(q).qp   *)
		    
		    W0 = (1/2) D( v.JvkT.Jvk.qp ,{q});
		    
		    W1 = (1/2) Transpose( D( (v.JvkT.Hat(Jok.qp) - v.JokT.Hat(Jvk.qp)).R0k  ,{q}) );
 			
 			W2 = (1/2) Transpose( D( v.JokT.R0k.EE.Transpose(R0k).Jok.qp ,{q}) );
 			
 			
 			% termini da d U   *)
 			%            -     *)
 			%            d q   *)
 			
 			Z0 = -JvkT.g0;
 			
 			%
 			if ( k == 1,
 			Z1 = -Transpose( { D( g0.R0k , {q}) } ),
 			Z1 = -Transpose(   D( g0.R0k , {q})   )
 			 );
 			 *)
 			 
 			Z1 = -Transpose(   D( g0.R0k , {q})   ); 
 			
 			Y0 = X0p - W0  + Z0;
 			Y1 = X1p - W1  + Z1;
 			Y2 = X2p - W2;
 			
 			Y  = StackCols(Y0, Y1, Y2);
 			Return(Y); 
 			 
 			 
 			 % test:   Return(X2p); *)	
	); 

Regressor(DHtable_, q_, qp_, v_, vp_, t_, g0_) 
	%Module({n, Yk, Y},
		
		n = Length(DHtable);
		Y = Regressor(DHtable, q, qp, v, vp, t, g0, 1);
		For( k = 2, k <= n, k++,
			
			Yk = Regressor(DHtable, q, qp, v, vp, t, g0, k);
			Y  = StackCols(Y, Yk);
		);
		
		Return(Y);
	);

Inertia(DHtable_, CGtable_, Masslist_, Tensortable_) 
	%Module(
		{Bmat,
			k, n = Length(DHtable),
		 R0k, Jk, Jvk, Jok, JvkT, JokT},
		 
		 Bmat = ZeroMatrix(n);
		
		For( k = 1, k <= n, k++,
			
			R0k  = RigidOrientation( DHFKine(DHtable, k) );
		    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
		    
		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
		    
		    Jvk =  Jk((1;;3, All));
		    Jok =  Jk((4;;6, All));
		    
		    JvkT = Transpose(Jvk);
		    JokT = Transpose(Jok);
		    
		    Bmat += Masslist((k)) JvkT.Jvk + JokT.R0k.Tensortable((k)).Transpose(R0k).Jok; 
		);
		Return(Bmat);
	);

InertiaBase(DHtable_, CGtable_, Masslist_, Tensortable_, {Tb0_, Tne_}) 
	%Module(
		{Bmat,
			k, n = Length(DHtable),
		 R0k, Jk, Jvk, Jok, JvkT, JokT},
		 
		 Bmat = ZeroMatrix(n);
		
		For( k = 1, k <= n, k++,
			
			R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
		    Jk  =  CGJacobBaseDyn(DHtable, CGtable, {Tb0, Tne}, k);
		    
		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
		    
		    Jvk =  Jk((1;;3, All));
		    Jok =  Jk((4;;6, All));
		    
		    JvkT = Transpose(Jvk);
		    JokT = Transpose(Jok);
		    
		    Bmat += Masslist((k)) JvkT.Jvk + JokT.R0k.Tensortable((k)).Transpose(R0k).Jok; 
		);
		Return(Bmat);
	);

InertiaToCoriolis(M_, q_, qp_)  
  %Module(
    {Cmat, i, j, k, n = Length(M)},

    % Brute force calculation *)
    Cmat = ZeroMatrix(n);

    For(i = 1, i <= n, ++i,
      For(j = 1, j <= n, ++j,
        For(k = 1, k <= n, ++k,
          Cmat((i,j)) += (1/2) * qp((k)) * (D(M((i,j)), q((k))) + D(M((i,k)), q((j))) - D(M((j,k)), q((i))))
	       )
         )
       );
   Cmat
  );

InertiaToCoriolisNotChristoffel(M_, q_, qp_)  
	%Module(
		{n = Length(M), Cmat, i, j, k},
		Cmat = ZeroMatrix(n);
		
	   	   
	   	   For(i = 1, i <= n, ++i,
             For(j = 1, j <= n, ++j,
                Cmat((i,j)) = Sum(
           	                  ( D( M((i, j)), q((k)) ) - (1/2) D( M((j, k)), q((i)) ) ) qp((k)), {k, 1, n}
           	                     )
             )
	   	   );
           Return(Cmat);	                  
	);
	
Gravitational(DHtable_, CGtable_, Masslist_, g0_) 
	%Module({Gvec, k, n = Length(DHtable)},
		
		Gvec = (ZeroMatrix(1, n))((1));
		
		For( k = 1, k <= n, k++,
			
			R0k  = RigidOrientation( DHFKine(DHtable, k) );
		    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
		    
		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
		    
		    Jvk =  Jk((1;;3, All));
		    % Jok =  Jk((4;;6, All)); *)
		    
		    JvkT = Transpose(Jvk);
		    % JokT = Transpose(Jok);  *)
		    
		    Gvec += -Masslist((k)) JvkT.g0; 
		);
		Return(Gvec);
		
	);

GravitationalBase(DHtable_, CGtable_, Masslist_, g0_, {Tb0_, Tne_}) 
	%Module({Gvec, k, n = Length(DHtable)},
		
		Gvec = (ZeroMatrix(1, n))((1));
		
		For( k = 1, k <= n, k++,
			
			R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
		    Jk  =  CGJacobBaseDyn(DHtable, CGtable, {Tb0, Tne}, k);
		    
		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
		    
		    Jvk =  Jk((1;;3, All));
		    % Jok =  Jk((4;;6, All)); *)
		    
		    JvkT = Transpose(Jvk);
		    % JokT = Transpose(Jok);  *)
		    
		    Gvec += -Masslist((k)) JvkT.g0; 
		);
		Return(Gvec);
		
	);
	
DynamicEquations(DHtable_, CGtable_, Masslist_, Tensortable_, g0_, q_, qp_, v_, vp_) 
	%Module({B, C, G, eqs},
		
		B = Inertia(DHtable, CGtable, Masslist, Tensortable);
		C = InertiaToCoriolis(B, q, qp);
		G = Gravitational(DHtable, CGtable, Masslist, g0);
		
		eqs = B.vp + C.v + G;
		
		Return(eqs);	
	);
	
DynamicEquationsBase(DHtable_, CGtable_, Masslist_, Tensortable_, gb_, q_, qp_, v_, vp_, {Tb0_, Tne_}) 
	%Module({B, C, G, eqs},
		
		B = InertiaBase(DHtable, CGtable, Masslist, Tensortable, {Tb0, Tne});
		C = InertiaToCoriolis(B, q, qp);
		G = GravitationalBase(DHtable, CGtable, Masslist, gb, {Tb0, Tne});
		
		eqs = B.vp + C.v + G;
		
		Return(eqs);	
	);
	
% parameters reshaping *)
TakeMoments(A_List) 
	%Module({v},
		v = {A((1,1)), -A((1,2)), -A((1,3)), A((2,2)), -A((2,3)), A((3,3))};
		
		Return(v);
		
	);

ExtractParameters(CGtable_, Masslist_, Tensortable_) 
	%Module({p, n, k},
		p = {};
		n = Length(CGtable);
		
		For( k = 1, k <= n, k++,
			
             p0k = {Masslist((k))};
             
             p1k = Masslist((k)) CGtable((k));
             
             p2k = TakeMoments( Tensortable((k)) - Masslist((k)) Hat(CGtable((k))).Hat(CGtable((k))) );
             				
			 pk = Join( p0k, p1k, p2k );
			 
			 p = Join(p, pk); 
		); 
		
		Return(p);
		
		
	);		  

SwitchAt(v1_List, v2_List, k_)  
	%Module({res,index},
		res = v1;
		index = k;
		res((index)) = v2((index));
		Return(res);
	);
	
% numerical derivative of a matrix w.r.t. q at numerical value q0 *)	

NDMatrix(S_, q_, q0_)  
	%Module({r, i, m, n},
		
		r = Length(q);
		{m, n} = Dimensions( S(q0) );
		dS = Table(ZeroMatrix(m,n),{r});
		
		
		For(i=1, i<=r, i++,
			
			% dS((i)) = ND( S /. Thread(Rule(q, SwitchAt(q0, q, i) )), q((i)), q0((i)) ); *)
			
			dS((i)) = ND( S(SwitchAt(q0, q, i)), q((i)), q0((i)) );
			
		);
		
		Return(dS);
	);

% the whole expression A(q1,q2,...) must be supplied *)
	
NullBasis(A_, independent_)  
	%Module({m, n, eye, dependent, depsel, indepsel, LHSmat, RHSvec, i},
		  
		   {m, n} = Dimensions(A);
		      eye = IdentityMatrix(n - m);
		dependent = Complement(Range(n), independent);
		 depsel   = SelectionMatrixColumns(n, dependent);
		 indepsel = SelectionMatrixColumns(n, independent);
		   LHSmat =  A.depsel;
		   RHSvec = -A.indepsel;
		      res = Inverse(LHSmat).RHSvec;
		      
		 For(i=1, i<=Length(independent), i++,
		 	
		 	  res = Insert(res, eye(( i )), {independent((i))} );
		 
		 );
		 
		 Return(res);
		      
		 
	);

% only the head of the matrix has to be supplied *)
	
NullBasisNumerical(A_, q0_, independent_)  
	%Module({m, n, Aq0},
		
		if ( MatrixQ(A(q0)) , Aq0 = A(q0), Aq0 = A@@q0);
		
		{m, n} = Dimensions( Aq0 );
		   eye = IdentityMatrix(n - m);
		   
		dependent = Complement(Range(n), independent);
		 depsel   = SelectionMatrixColumns(n, dependent);
		 indepsel = SelectionMatrixColumns(n, independent);
		   LHSmat =  Aq0.depsel;
		   RHSvec = -Aq0.indepsel;
		      res = PseudoInverse(LHSmat).RHSvec;
		      
		 For(i=1, i<=Length(independent), i++,
		 	
		 	  res = Insert(res, eye(( i )), {independent((i))} );
		 
		 );
		 
		 Return(res);
	);

Eye(n_)  IdentityMatrix(n);
	
% Constraint definitions *)

GraspMatrix(EP_) 
	%Module({G},
		G = BlockMatrix(
			{{ IdentityMatrix(3),     ZeroMatrix(3) },
            {           Hat(EP), IdentityMatrix(3) }}
		);
		
		Return(G);
	);

GlobalGraspMatrix(EPlist_) 
	%Module({G},
		G = StackCols@@GraspMatrix/@EPlist;
		
		Return(G);
	);
			
% RigidAdjoint(gab), where gab = RPToHomogeneous(Rab, pab), gives the desired twist trasformation for a change of pole given by pab and a change of orientation given by Rab*)
% see the twist notation functions *)

ConstraintMatrix(constrainttype_, axis_, Rbe_) 
	%Module({H, F, n,  p, RT, g, Adg},
		H = IdentityMatrix(6);
		
		p = {0, 0, 0};
		RT = Transpose(Rbe);
		
		g = RPToHomogeneous(RT, p);
		
		Adg = RigidAdjoint(g);
		
		if (constrainttype == "S", 
			H = Take(H, {1, 3}, {1, 6}); 
		);
		
		if (constrainttype == "R",  
		    F = Join( {0, 0, 0}, axis );
		    H = NullSpace({F});
		);
		
		if (constrainttype == "P",  
		    F = Join( axis, {0, 0, 0} );
		    H = NullSpace({F});
		);
		
		if (constrainttype == "C",  
		     H = IdentityMatrix(6);
		);
		
		if (constrainttype == "PC",  
		    H = Take(H, {3});
		);
		
		if (constrainttype == "PCWF",  
		    H = Take(H, {1, 3}, {1, 6});
		);
		
		if (constrainttype == "SF",  
		    H = H(({1, 2, 3, 6}, All));
		);
		
		Return(H.Adg);
	);

GlobalConstraintMatrix(list_) 
	%Module({consmats, H},
		
		consmats = (ConstraintMatrix @@ #)& /@ list;
		
		H        = BlockDiag(consmats);
		
		Return(H);
		
	);


FreeHandJacobian(list_) 
	%Module({Jacmats, fhJac},
		
		% each element of list must be {DHtable_i, vars_i, gp0i } *)
		
		Jacmats = DHJacobBase(#((1)) @@ #((2)), #((3)))& /@ list;
		fhJac   = BlockDiag(Jacmats);
		
		Return(fhJac);
		
	);
	
HandJacobian(fingerlist_, constraintlist_) 
	%Module({fhJac, H, hJac},
		
		fhJac = FreeHandJacobian(fingerlist);
		H     = GlobalConstraintMatrix(constraintlist);
		
		hJac  = H.fhJac;
		
		Return(hJac); 
	);	

ObjectJac( function_, pars_ ) 	
	%Module({res},
		
		res = BlockMatrix(
			{ 
			  { Eye(3)       , ZeroMatrix(3) },
			  { ZeroMatrix(3), function@@pars }
			}	
		);
		
		Return(res);
		
	);

CouplerJacobian( fstring_, pars_, constraintlist_, pointlist_ ) 
	%Module({H, G, GT, Jo,  fJac, fFKin, res },
		
		fFKin  = ToExpression(fstring <> "ToMat");
		fJac   = ToExpression(fstring <> "ToSpatialJac");
		
		H  = GlobalConstraintMatrix(constraintlist);
		
        G  = GlobalGraspMatrix( fFKin@@pars.(#) & /@ pointlist);
        GT = Transpose(G);
        
		Jo = ObjectJac(fJac, pars); 
 
	   res = H.GT.Jo;
	   
	   Return(res); 
		
	);

PfaffianMatrix(fingerlist_, fstring_, pars_, constraintlist_, pointlist_) 
	%Module({A11, A12, res},
		
		A11 =   HandJacobian(fingerlist, constraintlist);
		
		A12 = - CouplerJacobian(fstring, pars, constraintlist, pointlist);
		
		res = StackCols(A11, A12);
		
		Return(res);
	);

% Graphical Primitives *)

CreateGround({xmin_, xmax_, dx_}, {ymin_, ymax_, dy_}, colorplane_, colorframe_) 
	%Module({x, y, planepts, plane, e1s, e2s, e3s, Orig, frameS, scene},
		
		planepts = Table({x, y, 0}, {x, xmin, xmax, dx}, {y, ymin, ymax, dy});

        plane = Map( Graphics3D, {Map({colorplane, Line(#)} &, planepts), Map({colorplane, Line(#)} &, Transpose(planepts))});
        
        e1s = {5, 0, 0}; e2s = {0, 5, 0}; e3s = {0, 0, 5}; Orig = {0, 0, 0};
        
        frameS = Graphics3D({Thick, colorframe, Line({Orig, e1s}), Line({Orig, e2s}), Line({Orig, e3s})});

        scene = Append(plane, frameS);
        
        Return(scene);
		
	);

CreateFrame(g_, colorframe_) 
	%Module({R, p, x, y, z, scale, frame},
		scale = 5.0;
		
		R = RigidOrientation(g);
		p = RigidPosition(g);
		
		px = p + scale Flatten(TakeColumns(R, {1}));
		py = p + scale Flatten(TakeColumns(R, {2}));
		pz = p + scale Flatten(TakeColumns(R, {3}));
		
		frame = Graphics3D({Thick, colorframe, Line({p, px}), Line({p, py}), Line({p, pz})});
	
	    Return(frame);
	);
	
CreateJoint(g_, {height_, radius_}) 
	%Module({R, p, z, bottom, top, cylinder},
		
		R = RigidOrientation(g);
		p = RigidPosition(g);
		z = Flatten(TakeColumns(R, {3}));
		
		bottom = p - (height/2) z;
		top    = p + (height/2) z;
		
		cylinder = Graphics3D({
			                   Green, Opacity(0.25), Specularity(White, 30),
		                       Cylinder({bottom, top}, radius)
		                       }, Boxed -> False);
		
		Return(cylinder);	
	);
	
CreateLink(g_, {length_, radius_}) 
	%Module({R, p, x, start, end, cylinder},
		
		R = RigidOrientation(g);
		p = RigidPosition(g);
		x = Flatten(TakeColumns(R, {1}));
		
		start  = p - length x;
		end    = p;
		
		cylinder = Graphics3D({
			                   Orange, Opacity(0.25), Specularity(White, 30),
		                       Cylinder({start, end}, radius)
		                       }, Boxed -> False);
		
		Return(cylinder);	
	);
	
CreateRobot(DHtable_, {gpb_, gb0_}, q0_, a0_:20) 
	%Module({Tpb, Tb0, i, n, rlink, rjoint, hjoint, robot, res},
		Tpb = gpb;
		Tb0 = gb0;
		n   = Length(DHtable@@q0);
		rlink  = 2.5;
		rjoint = 5.0; 
		hjoint = 10.0;
		     % a0 = 20; *)
		
		% ruleq = Thread(Rule(q, q0)); *)
		
		link0  = {CreateLink(Tpb.Tb0.HomogeneousRotY(-Pi/2), {a0, rlink})};
		joint1 = {CreateJoint(Tpb, {hjoint, rjoint})};
		
		linklist = Table( 
			             CreateLink( 
			             	   DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)}, i 
			             	          ), 
			             		        { (DHtable@@q0)((i, 1)), rlink } 
			                       ),
			             {i, 1, n}
			             );
			             
		jointlist = Table( 
			             CreateJoint( 
			             	   DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)}, i-1 
			             	          ), 
			             		        { hjoint, rjoint } 
			                       ),
			             {i, 2, n}
			             );	             
		frameEE = {CreateFrame( DHFKine( DHtable@@q0, {Tpb.Tb0, Eye(4)} ), Magenta )};
			             
		res = Join(link0, joint1, linklist, jointlist, frameEE);	             
		
		Return(res);
	);

BBox(a_, b_, c_) 
	%Module({list},
		p1 = {a, b, -c, 1}; p2 = {-a, b, -c, 1}; p3 = {-a, -b, -c, 1}; p4 = {a, -b, -c, 1};
		p5 = {a, b,  c, 1}; p6 = {-a, b,  c, 1}; p7 = {-a, -b,  c, 1}; p8 = {a, -b,  c, 1};
	    
	    list1 = {p1, p2, p3, p4, p1}//Transpose;
	    list2 = {p5, p6, p7, p8, p5}//Transpose;
	    list3 = {p1, p5}//Transpose;
	    list4 = {p2, p6}//Transpose;
	    list5 = {p3, p7}//Transpose;
	    list6 = {p4, p8}//Transpose;
	    
	    list = {list1, list2, list3, list4, list5, list6};
	    
	    Return(list);
	
	);

EEllipsoid(a_, b_, c_) 
	%Module({nu, nv},
		
		nu = 20.; nv = 20.;
		
		x(u_, v_)  a Cos(v) Cos(u);
		y(u_, v_)  b Cos(v) Sin(u);
		z(u_, v_)  c Sin(v);
		
		r(u_, v_)  {x(u,v), y(u,v), z(u,v), 1};
		
		list1 = Table(r(u,v), {u, 0, 2. Pi, 2. Pi/nu}, {v, -Pi/2., Pi/2., Pi/nv});
		list2 = Transpose(list1);
		list1 = Transpose /@ list1;
		list2 = Transpose /@ list2;
		Return( Join(list1, list2) );
		
	);
	
CreateObject(g_, objectlist_, color_) 
	%Module({points, lines, res},
		
		points = Transpose(Drop(Dot(g,#),-1)) & /@ objectlist;
		
		lines = Line /@ points;
		
		res = Graphics3D({Thin, color, lines});
		
		Return(res);
		
	);
	
FindRedundantAnglesFromOrigin(origin_, DHtable_, {Tb0_, Tne_}, q_, qguess_) 
	%Module({g, p, f, L, gcons1, gcons2, gcons, cons, initguess, res},
		
		g = DHFKine(DHtable, {Tb0, Tne});
		p = RigidPosition(g);
		
		f      = {(1/2) q.q}; 
		% f      = {Max(q)};    *)
		
		gcons1  = origin - p;
		gcons2  = {};
		gcons   = Join(gcons1, gcons2);
		cons     = Thread(Equal(gcons, 0));
		
		L = Join(f, cons);
		
		initguess = Thread(List(q, qguess));
		
		res = q /. FindMinimum(L, initguess )((2)); 
		% res = q /. NMinimize(L, q )((2)); *)

        Return(res);		
	);