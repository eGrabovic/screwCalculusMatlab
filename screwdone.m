%% Select submatrices
%{
TakeRows[mat_?MatrixQ, part_] := Take[mat, part]

TakeColumns[mat_?MatrixQ, part_] := Take[mat, All, part]

TakeMatrix[mat_?MatrixQ, start:{startR_Integer, startC_Integer},
end:{endR_Integer, endC_Integer}] :=
	Take[mat, {startR, endR}, {startC, endC}] /;
	And @@ Thread[Dimensions[mat] >= start] && 
	And @@ Thread[Dimensions[mat] >= end]

SubMatrix[mat_List, start:{_Integer, _Integer}, dim:{_Integer,_Integer}] :=
	TakeMatrix[mat, start, start+dim-1]; 

ZeroMatrix[0, ___] := {};

ZeroMatrix[m_Integer, 0] := Table[{}, {m}];

ZeroMatrix[m_Integer,n_Integer] := 
	Normal[SparseArray[{},{m, n}]] /; m >= 0 && n>=0

ZeroMatrix[m_Integer] := ZeroMatrix[m, m] /; m >= 0 

SelectionMatrix[q_List, indexlist_List] :=
	Module[{m, eye},
		m = Length[q];
		eye = IdentityMatrix[m];
		res = eye[[indexlist]];
		Return[res];
	];

SelectionMatrix[m_Integer, indexlist_List] :=
	Module[{eye},
		eye = IdentityMatrix[m];
		res = eye[[indexlist]];
		Return[res];
	];

SelectionMatrixColumns[n_Integer, indexlist_List] :=	
	Module[{eye},
		eye = IdentityMatrix[n];
		res = Transpose[ eye[[indexlist]] ];
		Return[res];
	];

 % Extract the orientation portion from a homogeneous transformation *)
RigidOrientation(g_?MatrixQ) 
  %Module(
    {nr, nc},

    % Check to make sure that we were passed a square matrix *)
    if (Not(MatrixQ(g)) || ({nr, nc} = Dimensions(g); nr ~= nc) || nr < 3,
        Message(Screws::wrongDimensions, "RigidOrientation");
	Return Null;
    );

    % Extract the 3x3 upper left corner *)
    SubMatrix(g, {1,1}, {nc-1,nc-1})
  );

% Extract the orientation portion from a homogeneous transformation *)
RigidPosition(g_?MatrixQ) 
  %Module(
    {nr, nc},

    % Check to make sure that we were passed a square matrix *)
    if (Not(MatrixQ(g)) || ({nr, nc} = Dimensions(g); nr ~= nc) || nr < 3,
        Message(Screws::wrongDimensions, "RigidPosition");
	Return Null;
    );

    % Extract the upper left column *)
    Flatten(SubMatrix(g, {1, nc}, {nc - 1 ,1}))
  );
%}

%% Transformations
%{
(* Rotations about a single axis*)
RotX[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{1,  0, 0},
          {0, ca, -sa}, 
          {0, sa,  ca}
          }
  ];

RotY[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{ca,  0, sa},
          {0,   1,  0}, 
          {-sa, 0,  ca}
          }
  ];

RotZ[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{ca, -sa, 0},
          {sa,  ca, 0}, 
          {0,    0, 1}
          }
  ];
%}

%% Already existent
%{
ATan2[y_, x_] := ArcTan[x,y]; --> atan2()

Magnitude(v_)  --> norm()
%Module(
  {},
  
  if (Not(VectorQ(v)),
    Message(Screws::wrongDimensions, "Vector");
    Return(Null);
  );

  Sqrt(v.v)
);

function BlockDiag(list_)  --> diag()
	%Module({r, mlist, nlist, m, n, i},
		r = Length(list);
		mlist = Dimensions(#)((1))& /@ list;
		nlist = Dimensions(#)((2))& /@ list;
		m = Plus @@ mlist;
		n = Plus @@ nlist;
		res = ZeroMatrix(m, n);
		
		res(( 1;;mlist((1)), 1;;nlist((1)) ))= list((1));
		
		
		For( i=2, i <= r, ++i,
			res(( Plus@@Take(mlist, i-1)+1 ;; Plus@@Take(mlist, i), Plus@@Take(nlist, i-1)+1 ;; Plus@@Take(nlist, i) )) = list((i))
		);
		
		Return(res);
	);



Eye(n_)  IdentityMatrix(n); --> eye()
%}

%% Skews and anti
%{
% Generate a skew symmetric matrix from an axis*)

Hat(w_)   AxisToSkew(w);  % synonym *) 

Skew(w_)   AxisToSkew(w); % synonym *)

AxisToSkew(omega_?VectorQ) 
  %Module(
    {},
    % Check to make sure the dimensions are okay *)
    if (Not(VectorQ(omega)) || Length(omega) ~= 3,
      Message(Screws::wrongDimension);
      Return Null;
    );

    % Return the appropriate matrix *)
    {{ 0,          -omega((3)),  omega((2))},
     { omega((3)), 0,           -omega((1))},
     {-omega((2)), omega((1)),  0          }}
  );

% Generate an axis from a skew symmetric matrix *)
HatInv(S_)   SkewToAxis(S); % synonyms *)

UnSkew(S_)   SkewToAxis(S);  % synonym *) 

SkewToAxis(S_) 
  %Module(
    {},
    % First check to make sure we have a skew symmetric matrix *)
    if (Not(skewQ(S)) || Dimensions(S) ~= {3,3},
      Message(Screws::wrongDimension);
      Return Null
    );

    % Now extract the appropriate component *)
    {S((3,2)), S((1,3)), S((2,1))}
  );

SkewExp(S_?skewQ,theta_:1) 
  %Module(
    {n = Dimensions(S)((1))},

    % Use Rodrigues's formula *)
    IdentityMatrix(3) + Sin(theta) S + (1 - Cos(theta)) S.S
  );


%}
%% Checks
%{
  function skewQ(mat_)  
%   %Module({nr, nc, zmat},

    % First check to make sure that this is square matrix *)
    if Not(MatrixQ(mat)) || ({nr, nc} = Dimensions(mat); nr ~= nc)
	    Message(notSquare);    
        Return(False);
    end

    % Check to see if A = -A^T *)
    zmat = mat + Transpose(mat);
%     Return( And @@ Map(TrueQ(Simplify(#) == 0)&, Flatten(zmat)))
end

 function RotationQ(mat_)  
%   %Module({nr, nc, zmat},

    % First check to make sure that this is square matrix *)
    if Not(MatrixQ(mat)) || {nr, nc} == Dimensions(mat); nr ~= nc
        Message(notSquare);    
        Return(False);
    end

    % Check to see if R^T R = Identity *)
    zmat = Simplify(mat . Transpose(mat)) - IdentityMatrix(nr);
    %Return( And @@ Map(TrueQ(Simplify(#) == 0)&, Flatten(zmat)))
 end


%}


%% Quaternions
%{


QLength(curveOnSphere_List)  
  Sum(ArcCos(curveOnSphere((i)) . curveOnSphere((i+1))),
        {i,1,Length(curveOnSphere)-1})

QuatInv(q_List)   {q((1)), - q((2)), -q((3)), - q((4))}

QuatProd[q_List,p_List] :=
{q[[1]]*p[[1]] - q[[2]]*p[[2]] - q[[3]]*p[[3]] - q[[4]]*p[[4]],
 q[[1]]*p[[2]] + q[[2]]*p[[1]] + q[[3]]*p[[4]] - q[[4]]*p[[3]],
 q[[1]]*p[[3]] + q[[3]]*p[[1]] + q[[4]]*p[[2]] - q[[2]]*p[[4]],
 q[[1]]*p[[4]] + q[[4]]*p[[1]] + q[[2]]*p[[3]] - q[[3]]*p[[2]]} /;
   Length[q] == Length[p] == 4


normalize(lst_List)   %Module({norm = lst . lst,eps=10^(-14)},
if (NumberQ(norm),
if (norm<eps, lst, N(lst/Sqrt(norm)) ),
lst/Sqrt(norm)))  


proj4 (pt4D_,angle_:0)   Chop({pt4D((1)),pt4D((2)), 
                   Cos(angle)*pt4D((3)) + Sin(angle)*pt4D((4))})//N

 % QuatToMat(q_List)  
%Module({q0= q((1)), q1=q((2)), q2 = q((3)), q3 = q((4))},
   %Module({d23 = 2 q2 q3, d1 = 2 q0 q1,
          d31 = 2 q3 q1, d2 = 2 q0 q2,
          d12 = 2 q1 q2, d3 = 2 q0 q3,
          q0sq = q0^2, q1sq = q1^2, q2sq = q2^2, q3sq = q3^2},
          {{q0sq + q1sq - q2sq - q3sq, d12 - d3, d31 + d2},
           {d12 + d3, q0sq - q1sq + q2sq - q3sq, d23 - d1},
           {d31 - d2, d23 + d1, q0sq - q1sq - q2sq + q3sq}} )); *)
          
QuatToMat(q_List)  
%Module( {beta0 = q((1)), betas = q(({2,3,4})),
         Im, 
         hatbetas,
         mat},
		 Im = IdentityMatrix(3);
		 hatbetas = Skew(betas);
         mat = Im + 2  hatbetas . ( beta0 Im + hatbetas);
         Return(mat);
);

MatToQuat(mat_List)  
   %Module({q0,q1,q2,q3,trace,s,t1,t2,t3},
         trace = Sum(mat((i,i)),{i,1,3});
         if (trace > 0,
             s = Sqrt(trace + 1.0);
               q0 = s/2;
               s = 1/(2 s);
               q1 = (mat((3,2)) - mat((2,3)))*s;
               q2 = (mat((1,3)) - mat((3,1)))*s;
               q3 = (mat((2,1)) - mat((1,2)))*s,
            if ((mat((1,1)) >= mat((2,2))) &&
                (mat((1,1)) >= mat((3,3))), % i=0,  j = 1, k = 2 *)
                  s = Sqrt(mat((1,1)) - mat((2,2)) - mat((3,3)) + 1.0);
                  q1 = s/2; s = 1/(2 s);
                  q0 = (mat((3,2)) - mat((2,3)))*s;
                  q2= (mat((2,1)) + mat((1,2)))*s;
                  q3 = (mat((1,3)) + mat((3,1)))*s,
               if ((mat((1,1)) < mat((2,2))) &&
                (mat((1,1)) >= mat((3,3))), % i=1,  j = 2, k = 0 *)
                  s = Sqrt(mat((2,2)) - mat((3,3)) - mat((1,1)) + 1.0);
                  q2 = s/2; s = 1/(2 s);
                  q0 = (mat((1,3)) - mat((3,1)))*s;
                  q3= (mat((3,2)) + mat((2,3)))*s;
                  q1 = (mat((2,1)) + mat((1,2)))*s,
             % Else: (mat((1,1)) < mat((2,2))) && (mat((1,1)) < mat((3,3)))  *)
             % i=2,  j = 0, k = 1 *)
            s = Sqrt(mat((3,3)) - mat((1,1)) - mat((2,2)) + 1.0);
                  q3 = s/2; s = 1/(2 s);
                  q0 = (mat((2,1)) - mat((1,2)))*s;
                  q1= (mat((1,3)) + mat((3,1)))*s;
                  q2 = (mat((3,2)) + mat((2,3)))*s)));
                  normalize(N({q0,q1,q2,q3})))

QuatVectPart(mat_List)   Take( MatToQuat(mat), -3);

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
 			


%}

%% Parametrizations of SO(3) - Forward Map *)
%{

EulZXZToMat(alpha_, beta_, gamma_) 
	%Module( {},
	RotZ(alpha).RotX(beta).RotZ(gamma)
	
	);

EulZYZToMat(alpha_, beta_, gamma_) 
	%Module( {},
	RotZ(alpha).RotY(beta).RotZ(gamma)
	
	);

EulZYXToMat(alpha_, beta_, gamma_) 
	%Module( {},
	RotZ(alpha).RotY(beta).RotX(gamma)
	
	);
	
RPYToMat(alpha_, beta_, gamma_)  EulZYXToMat(alpha, beta, gamma);

% with theta \in (0, Pi) *)
MatToEulZYZ(R_)  
	%Module({phi, theta, psi, thetaiszero, thetaisPi},
	        
	% Check the singularity of representation ZYZ *)        
	thetaiszero = Abs( R((3,3)) - 1) < 10^(-10);
	thetaisPi = Abs(R((3,3)) + 1) < 10^(-10);     
	        
	% In cases of singularity we set arbitrarily psi = 0 *)        
	        
	if ( thetaiszero,  
	     phi = ATan2( R((2,1)) , R((1,1)) );
	     theta = 0;
	     psi = 0;
	     );
	  
	if ( thetaisPi,
	     phi = ATan2( R((2,1)), R((1,1)) );
	     theta = Pi;
	     psi = 0;
	     );  
	   
	if ( !(thetaiszero || thetaisPi),
	     phi = ATan2(R((2,3)), R((1,3)));
	     theta = ATan2(Sqrt( R((1,3))^2 + R((2,3))^2 ), R((3,3)) );
	     psi = ATan2(R((3,2)), -R((3,1)));
	     );
	     
	     Return({phi, theta, psi});
	     
	);

(* with theta \in (-Pi/2, Pi/2) *)
MatToEulZYX[R_] :=
	Module[{phi, theta, psi, thetaisplushalfPi, thetaisminushalfPi},
	        
	(* Check the singularity of representation ZYX *)        
	thetaisplushalfPi = Abs[ R[[3,1]] + 1] < 10^(-10);
	thetaisminushalfPi = Abs[R[[3,1]] - 1] < 10^(-10);     
	        
	(* In cases of singularity we set arbitrarily psi = 0 *)        
	        
	If[ thetaisplushalfPi,  
	     phi = ATan2[ R[[2,3]] , R[[1,3]] ];
	     theta = Pi/2;
	     psi = 0;
	     ];
	  
	If[ thetaisminushalfPi,
	     phi = ATan2[ -R[[2,3]], -R[[1,3]] ];
	     theta = -Pi/2;
	     psi = 0;
	     ];  
	   
	If[ !(thetaisplushalfPi || thetaisminushalfPi),
	     phi = ATan2[R[[2,1]], R[[1,1]]];
	     theta = ATan2[ -R[[3,1]], Sqrt[ R[[3,2]]^2 + R[[3,3]]^2 ] ];
	     psi = ATan2[ R[[3,2]], R[[3,3]] ];
	     ];
	     
	     Return[{phi, theta, psi}];
	     
	];(* with theta \in (-Pi/2, Pi/2) *)
MatToEulZYX[R_] :=
	Module[{phi, theta, psi, thetaisplushalfPi, thetaisminushalfPi},
	        
	(* Check the singularity of representation ZYX *)        
	thetaisplushalfPi = Abs[ R[[3,1]] + 1] < 10^(-10);
	thetaisminushalfPi = Abs[R[[3,1]] - 1] < 10^(-10);     
	        
	(* In cases of singularity we set arbitrarily psi = 0 *)        
	        
	If[ thetaisplushalfPi,  
	     phi = ATan2[ R[[2,3]] , R[[1,3]] ];
	     theta = Pi/2;
	     psi = 0;
	     ];
	  
	If[ thetaisminushalfPi,
	     phi = ATan2[ -R[[2,3]], -R[[1,3]] ];
	     theta = -Pi/2;
	     psi = 0;
	     ];  
	   
	If[ !(thetaisplushalfPi || thetaisminushalfPi),
	     phi = ATan2[R[[2,1]], R[[1,1]]];
	     theta = ATan2[ -R[[3,1]], Sqrt[ R[[3,2]]^2 + R[[3,3]]^2 ] ];
	     psi = ATan2[ R[[3,2]], R[[3,3]] ];
	     ];
	     
	     Return[{phi, theta, psi}];
	     
	];

MatToRPY(R_)   MatToEulZYX(R);	


%}

%% Utils
%{
FramesToVect(list_)  
	%Module({ R  = list((1)),
			 Rd = list((2)),
			 hn,  hs,    ha,
			  nd,  sd,    ad,
			 eO },
			
			 hn = Hat(R((All, 1)));    
			 hs = Hat(R((All, 2)));   
			 ha = Hat(R((All, 3)));
			
		   	 nd = Rd((All, 1)); 
		   	 sd = Rd((All, 2)); 
		   	 ad = Rd((All, 3));
			
			 eO = (1/2) (hn.nd + hs.sd + ha.ad)
	);

MatToVect(R_)   
	%Module({axis, theta},
	{axis, theta} = RotationParam(R);
	axis*Sin(theta)
	)/; Length(R)==3  

MatToVect(list_)	 
	%Module({ R  = list((1)),
			 Rd = list((2)),
			 hn,  hs,    ha,
			  nd,  sd,    ad,
			 eO },
			
			 hn = Hat(R((All, 1)));    
			 hs = Hat(R((All, 2)));   
			 ha = Hat(R((All, 3)));
			
		   	 nd = Rd((All, 1)); 
		   	 sd = Rd((All, 2)); 
		   	 ad = Rd((All, 3));
			
			 eO = (1/2) (hn.nd + hs.sd + ha.ad)
	)/; Length(list)==2

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
	
%}

%% Rot Param
%{
RotationParam(R_)  
  %Module(
    {nr, nc},

    % Check to make sure that our input makes sense *)
    if (Not(MatrixQ(R)) || ({nr, nc} = Dimensions(R); nr ~= nc) || nr ~= 3,
        Message(Screws::wrongDimensions, "RotationAxis");
	Return(Null);
    );
     
     
    t = (Sum(R((i,i)),{i, nr})-1)/2;
    if (t<-1, t=-1;);
    if (t>1,t=1;);

    theta = ArcCos(t);

    if (theta ~= 0,
       axis=RotationAxis(R, theta);,
       axis = Table(0, {nc});
       theta=0;
    );

    Return({axis, theta}); 
 );

% Find the axis of a rotation matrix *)
RotationAxis(R_, theta_)  
  %Module(
    {nr, nc},

    % Check to make sure that our input makes sense *)
    if (Not(MatrixQ(R)) || ({nr, nc} = Dimensions(R); nr ~= nc),
        Message(Screws::wrongDimensions, "RotationAxis");
	Return(Null);
    );

    if (theta<0 || theta>Pi,
        Message(Screws::wrongTheta, "RotationAxis");
	Return(Null);
    );
 
    if (theta==Pi,
      axis=NullSpace(R-IdentityMatrix(3))((1));
      axis=axis/Magnitude(axis);
    ,
      axis={R((3,2))-R((2,3)),R((1,3))-R((3,1)),R((2,1))-R((1,2))}/(2*Sin(theta));
    );
    Return(axis);
);

%}

%% Homogeneous representation
%{

HomogeneousToVector(p_)  
	Block({},
		Take(p, 3)
	);

% Convert a point into homogeneous coordinates *)
PointToHomogeneous(p_)  
  Block({},
    % Check to make sure the dimensions of the args make sense *)
    % if (Not(VectorQ(p)), Message(Screws::notVector, "PointToHomogeneous")); *)

    % Now put everything together into a homogeneous vector *)
    Append(p, 1)
  ); 

% Convert a vector into homogeneous coordinates *)
VectorToHomogeneous(p_)  
  Block({},
    % Check to make sure the dimensions of the args make sense *)
    % if (Not(VectorQ(p)), Message(Screws::notVector, "VectorToHomogeneous")); *)

    % Now put everything together into a homogeneous vector *)
    Append(p, 0)
  );

% Convert a rotation + a translation to a homogeneous matrix *)
RPToHomogeneous(R_, p_)  
  %Module(
    {n},

    % Check to make sure the dimensions of the args make sense *)
    if (Not(VectorQ(p)) || Not(MatrixQ(R)) ||
       (n = Length(p); Dimensions(R) ~= {n, n}),
	Message(Screws::wrongDimensions, "RPToHomogeneous:");
    );

    % Now put everything together into a homogeneous transformation *)
    
    BlockMatrix({{R,       Transpose({p})},
                 {ZeroMatrix(1,3),  {{1}}}
                 }
    )
    
  );  

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
%}

%% Rigid
%{
% Calculate the inverse rigid transformation *)

RigidInverse(g_?MatrixQ)   
  %Module(
    {R = RigidOrientation(g), p = RigidPosition(g)},
    RPToHomogeneous(Transpose(R), -Transpose(R).p)
  );  

 % Figure out the dimension of a twist (private) *)
xidim(xi_?VectorQ)  
  %Module(
    {l = Length(xi), n},

    % Check the dimensions of the vector to make sure everything is okay *)
    n = (Sqrt(1 + 8l) - 1)/2;
    if (Not(IntegerQ(n)),
      Message(Screws::wrongDimensions, "xidim");
      Return 0;
    );
    n
);

% Extract the linear portion of a twist (private) *)
xitov(xi_?VectorQ)  
  %Module(
    {n = xidim(xi)},

    % Make sure that the vector had a reasonable length *)   
    if (n == 0, Return Null);

    % Extract the linear portion of the twist *)
    % SetPrecision(Take(xi, n),PRECISION) *)
    Take(xi, n)
  );




%}

%% done by grabovic
%{
% Adjoint matrix calculation *)
RigidAdjoint(g_?MatrixQ)   -> adjoint
  %Module(
    {R = RigidOrientation(g), p = RigidPosition(g)},
    
    BlockMatrix({{R,       AxisToSkew(p) . R},
                 {ZeroMatrix(3,3),  R}
                 }
    )
    
  );
  
% Inverse adjoint matrix calculation *)
InverseRigidAdjoint(g_?MatrixQ)   -> adointInv
  %Module(
    {RT = Transpose(RigidOrientation(g)), p = RigidPosition(g)},
    
    BlockMatrix({{RT,       -RT.AxisToSkew(p)},
                 {ZeroMatrix(3,3),  RT}
                 }
    )
    
  );  

% Construct the Body Jacobian for a robot with any no. of links *)			
BodyJacobian(args__, gst0_)  -> BodyJac 
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

% Construct the Spatial Jacobian for a robot with any no. of links *)
SpatialJacobian(args__, gst0_)   -> spatialJac
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

% Gives the homogeneous matrix *)
ForwardKinematics(args__, gst0_)  -> FWKin
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
%}

%% To be updated or done
%{
% Check to see if a matrix is a twist matrix *)
%! Not implemented !*)
TwistMatrixQ(A_)   MatrixQ(A);

% Convert a homogeneous matrix to a twist *)
%! This only works in dimensions 2 and 3 for now !*)
HomogeneousToTwist(A_)  
  %Module(
    {nr, nc},

    % Check to make sure that our input makes sense *)
    if (Not(MatrixQ(A)) || ({nr, nc} = Dimensions(A); nr ~= nc),
        Message(Screws::wrongDimensions, "HomogeneousToTwist");
	Return Null;
    );

    % Make sure that we have a twist and not a rigid motion *)
    if (A((nr,nc)) ~= 0,
        Message(Screws::notTwistMatrix, "HomogeneousToTwist");
	Return Null;
    );

    % Extract the skew part and the vector part and make a vector *)
    Join(
      Flatten(SubMatrix(A, {1, nr}, {nr - 1, 1})),
      SkewToAxis( SubMatrix(A, {1, 1}, {nr - 1 ,nr - 1}) )
    )
  );



% Convert a twist to homogeneous coordinates *)
TwistToHomogeneous(xi_?VectorQ)  
  %Module(
    {w = xitow(xi), v = xitov(xi), R, p},
    
    % Make sure that we got a real twist *)
    if (w == Null || v == NULL, Return Null);

    % Now put everything together into a homogeneous transformation *)
    BlockMatrix({{AxisToSkew(w),   Transpose({v})},
                 {ZeroMatrix(1,3), {{0}} }
                 }
    )
  );  


%}

%% Twists
%{
% Calculation of the error twist xi_err and error angle theta_err such that:
% TwistExp(xi_err, theta_err).gStart = gEnd
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


 % Calculation of the error twist xi_err and error angle theta_err such that:
%     gStart.TwistExp(xi_err, theta_err) = gEnd
    
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

% Find the magnitude of a twist *)
TwistMagnitude(xi_?VectorQ)   
  %Module(
    {v, w},
    {v, w} = Partition(xi, 3);
    if ((MatchQ(w,{0,0,0}) || MatchQ(w, {{0},{0},{0}})), 
      Sqrt(v.v), Sqrt(w.w))
  );

% Find the twist which generates a rigid motion - OLD VERSION *)
% RigidTwist(g_?MatrixQ)  
   %Module(
    {R, p, axis, v, theta, w},

    % Make sure the dimensions are okay *)
    %! Missing !*)

    % Extract the appropriate pieces of the homogeneous transformation *)
    R = RigidOrientation(g);
    p = RigidPosition(g);

    % Now find the axis from the rotation *)    
    %w = RotationAxis(R);
    *theta = RotationAngle(R);*)

    %Lagt till ist\(ADoubleDot)llet f\(ODoubleDot)r ovan*)
    {w,theta}=RotationParam(R);

    % Split into cases depending on whether theta is zero *)
    if (theta == 0,
      theta = Magnitude(p);
      if (theta == 0,  
        Return({{0,0,0,0,0,0},0});
        );
      v = p/theta;,
    % else *)
      % Solve a linear equation to figure out what v is *)   
      v = LinearSolve(
        (IdentityMatrix(3)-Outer(Times,w,w)) Sin(theta) +
        Skew(w) (1 - Cos(theta)) + Outer(Times,w,w) theta,
      p);
    );
	Return({Flatten({v, w}),theta});
 
  );
*)

(* Find the twist which generates a rigid motion - NEW VERSION *)
RigidTwist[g_?MatrixQ] :=
  Module[
    {R, p, axis, v, theta, w, Ainv},

    (* Make sure the dimensions are okay *)
    (*! Missing !*)

    (* Extract the appropriate pieces of the homogeneous transformation *)
    R = RigidOrientation[g];
    p = RigidPosition[g];

    (* Now find the axis from the rotation *)    
    (*w = RotationAxis[R];
    *theta = RotationAngle[R];*)

    {w,theta}=RotationParam[R];
    hatw = Hat[w];
     
    (* Split into cases depending on whether theta is zero *)
    If[theta == 0,
      theta = Magnitude[p];
      If[theta == 0,  
        Return[{{0,0,0,0,0,0},0}];
        ];
      v = p/theta;,
    (* else *)
      (* Solve a linear equation to figure out what v is *)   
      Ainv = IdentityMatrix[3]/theta - (1/2) hatw + (1/theta - (1/2) Cot[theta/2]) MatrixPower[hatw, 2];
      
      v = Ainv.p;
    ];
	Return[{Flatten[{v, w}],theta}];
 
  ];

% TwistExp(xi_?VectorQ, theta_:1)  
  %Module(
    {w = xitow(xi), v = xitov(xi), R, p},
      
    % Make sure that we got a real twist *)
    if (w == Null || v == NULL, Return Null);

    % Use the exponential formula from MLS *)
    if  ((MatchQ(w,{0,0,0}) || MatchQ(w, {{0},{0},{0}})),
      R = IdentityMatrix(3);
      p = v * theta;,
     % else *)
      ws=Skew(w);
      R = SkewExp(ws, theta);
      p = (IdentityMatrix(3) - R) . (ws . v) + w (w.v) theta;
    );
    RPToHomogeneous(R, p)
  );

  (* Build a twist from a screw *)
ScrewToTwist[Infinity, q_, w_] := Join[w, {0,0,0}];

ScrewToTwist(h_, q_, w_)   Join(-AxisToSkew(w) . q + h w, w)


% Find the pitch associated with a twist in R^3 *)
TwistPitch(xi_?VectorQ)   
  %Module(
    {v, w},
    {v, w} = Partition(xi, 3);
    v . w / w.w
  );

WrenchPitch(xi_?VectorQ)   Null;

% Find the axis of a twist *)
TwistAxis(xi_?VectorQ)   
  %Module(
    {v, w},
    {v, w} = Partition(xi, 3);
    if ((MatchQ(w,{0,0,0}) || MatchQ(w, {{0},{0},{0}})), 
     {0, v / Sqrt(v.v)}, {AxisToSkew(w) . v / w.w, (w / w.w)})
  );

 % Gives Xi 6 vector given a point on axis and axis unit vector for a Revolute Joint *)
RevoluteTwist(q_,w_)  Flatten({Cross(q,w),w});


% Gives Xi 6 vector given a point on axis and axis unit vector for a Prismatic Joint *)
PrismaticTwist(q_,w_)  Flatten({w, {0,0,0}});


%}

%% Reverse order?
%{
MakeQuat[n_List:{0,0,1},angle_:0] :=
     Module[{c = Cos[angle/2], s = Sin[angle/2]},
     {c, n[[1]]*s, n[[2]]*s, n[[3]]*s}//N]

MakeRot(n_List:{0,0,1},angle_:0)  
%Module({c = Cos(angle)//N, s = Sin(angle)//N, cm},
cm = 1 - c;
{{c + cm*n((1))^2, cm*n((2))*n((1)) - s*n((3)), cm*n((3))*n((1)) + s*n((2))},
{cm*n((1))*n((2)) + s*n((3)), c + cm*n((2))^2,  cm*n((3))*n((2)) - s*n((1))},
{ cm*n((1))*n((3)) - s*n((2)), 
              cm*n((2))*n((3)) + s*n((1)), c + cm*n((3))^2}}//N)
% ::Subsection:: *)
%Warning: these are the reverse of conventional argument order:*)
%change order once all is consistent to (angle, nhat).*)

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

%}
%% DH
%{
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
%% Jacob
%{

DHJacob0Dyn(DHtable_, k_) 	
% 	%Module({ i, n = Length(DHtable), Hk,
% 		     type, z, r, 
%              jv, jo, j, Jk, J },
% 		     
% 		Jk  = {};     
% 		Hk = DHFKine(DHtable, k);
% 		For(i = 1, i <= k, i++,
% 			
% 			type = DHtable((i, 5));
% 			   z = (RigidOrientation(DHFKine(DHtable, i-1)))((All, 3));
% 			   r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, i-1)));
% 			
% 			if ( type == "P",
% 				
% 				% Prismatic Joint *)
% 				jv = z;
% 				jo = {0, 0, 0};
% 				j  = Join(jv, jo) ; ,
% 			
% 			    % Revolute Joint *)
% 			    jv = Hat(z).r;
% 			    jo = z;
% 			    j  = Join(jv, jo); 
% 		      );
% 		Jk = Append(Jk, j);
% 		
% 		);
% 		
% 	    J = StackCols( Transpose(Jk), ZeroMatrix(6, n-k) );
% 	    Return(J);
% 	);

% DHJacobBaseDyn(DHtable_, {Tb0_, Tne_}, k_) 	
% 	%Module({ i, n = Length(DHtable), Hk,
% 		     type, z, r, 
%              jv, jo, j, Jk, J },
% 		     
% 		Jk  = {};    
% 		Hk  = DHFKine(DHtable, {Tb0, Tne}, k);
% 	
% 		For(i = 1, i <= k, i++,
% 			
% 			type = DHtable((i, 5));
% 			   z = (RigidOrientation(DHFKine(DHtable, {Tb0, Tne}, i-1)))((All, 3));
% 			   r = RigidPosition(Hk) - (RigidPosition(DHFKine(DHtable, {Tb0, Tne}, i-1)));
% 			
% 			if ( type == "P",
% 				
% 				% Prismatic Joint *)
% 				jv = z;
% 				jo = {0, 0, 0};
% 				j  = Join(jv, jo) ; ,
% 			
% 			    % Revolute Joint *)
% 			    jv = Hat(z).r;
% 			    jo = z;
% 			    j  = Join(jv, jo); 
% 		      );
% 		Jk = Append(Jk, j);
% 		
% 		);
% 		
% 	    J = StackCols( Transpose(Jk), ZeroMatrix(6, n-k) );
% 	    Return(J);
% 	);

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

%}
%% Wrench
%{
WrenchAxis(xi_?VectorQ)   Null;


WrenchMagnitude(xi_?VectorQ)   Null;
%}

%% Other
%{
(* numerical derivative of a matrix w.r.t. q at numerical value q0 *)	

NDMatrix[S_, q_, q0_] :=
	Module[{r, i, m, n},
		
		r = Length[q];
		{m, n} = Dimensions[ S[q0] ];
		dS = Table[ZeroMatrix[m,n],{r}];
		
		
		For[i=1, i<=r, i++,
			
			(* dS[[i]] = ND[ S /. Thread[Rule[q, SwitchAt[q0, q, i] ]], q[[i]], q0[[i]] ]; *)
			
			dS[[i]] = ND[ S[SwitchAt[q0, q, i]], q[[i]], q0[[i]] ];
			
		];
		
		Return[dS];
	];


ForceCloseQFrames[qfrms_List] :=
  Module[{lastfrm = qfrms[[1]], thisfrm},
    Table[thisfrm = qfrms[[i]];
          lastfrm =
             If[thisfrm . lastfrm >= 0, thisfrm, - thisfrm],
             {i,1,Length[qfrms]}]]

ForceClose2DQFrames(qrows_List)  
  %Module({lastrow = qrows((1)), thisrow},
    Table(thisrow = ForceCloseQFrames(qrows((i)));
          lastrow = 
             if (thisrow((1)) . lastrow((1)) >= 0, thisrow, - thisrow),
              {i,1,Length(qrows)}))
              
ForceClose3DQFrames(qplanes_List) %Module({lastplane=qplanes((1)),thisplane},Table(thisplane=ForceClose2DQFrames(qplanes((i)));
lastplane=if (thisplane((1)).lastplane((1))>=0,thisplane,-thisplane),{i,1,Length(qplanes)}))

%}
%% CG
%{

CGJacobBaseDyn(DHtable_, CGtable_, {Tb0_, Tne_}, k_) 
% 	%Module({R0k, pkck, M, DHJacob, CGJacob},
% 		R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
% 		pkck = CGtable((k));
% 		M = BlockMatrix( {
% 			               {   IdentityMatrix(3),     -Hat(R0k.pkck) },
% 						   {     ZeroMatrix(3,3),  IdentityMatrix(3) }
% 						  }
% 			           );
% 		DHJacob = DHJacobBaseDyn(DHtable, {Tb0, Tne}, k);
% 		CGJacob = M.DHJacob;
% 		Return(CGJacob);	           
% 	);

CGJacob0Dyn(DHtable_, CGtable_, k_) 
% 	%Module({R0k, pkck, M, DHJacob, CGJacob},
% 		R0k  = RigidOrientation( DHFKine(DHtable, k) );
% 		pkck = CGtable((k));
% 		M = BlockMatrix( {
% 			               {   IdentityMatrix(3),     -Hat(R0k.pkck) },
% 						   {     ZeroMatrix(3,3),  IdentityMatrix(3) }
% 						  }
% 			           );
% 		DHJacob = DHJacob0Dyn(DHtable, k);
% 		CGJacob = M.DHJacob;
% 		Return(CGJacob);	           
% 	);

%}

%% Constraints
%{

 PfaffianMatrix(fingerlist_, fstring_, pars_, constraintlist_, pointlist_) 
% 	%Module({A11, A12, res},
% 		
% 		A11 =   HandJacobian(fingerlist, constraintlist);
% 		
% 		A12 = - CouplerJacobian(fstring, pars, constraintlist, pointlist);
% 		
% 		res = StackCols(A11, A12);
% 		
% 		Return(res);
% 	);

CouplerJacobian( fstring_, pars_, constraintlist_, pointlist_ ) 
% 	%Module({H, G, GT, Jo,  fJac, fFKin, res },
% 		
% 		fFKin  = ToExpression(fstring <> "ToMat");
% 		fJac   = ToExpression(fstring <> "ToSpatialJac");
% 		
% 		H  = GlobalConstraintMatrix(constraintlist);
% 		
%         G  = GlobalGraspMatrix( fFKin@@pars.(#) & /@ pointlist);
%         GT = Transpose(G);
%         
% 		Jo = ObjectJac(fJac, pars); 
%  
% 	   res = H.GT.Jo;
% 	   
% 	   Return(res); 
% 		
% 	);

ObjectJac( function_, pars_ ) 	
% 	%Module({res},
% 		
% 		res = BlockMatrix(
% 			{ 
% 			  { Eye(3)       , ZeroMatrix(3) },
% 			  { ZeroMatrix(3), function@@pars }
% 			}	
% 		);
% 		
% 		Return(res);
% 		
% 	);

HandJacobian(fingerlist_, constraintlist_) 
% 	%Module({fhJac, H, hJac},
% 		
% 		fhJac = FreeHandJacobian(fingerlist);
% 		H     = GlobalConstraintMatrix(constraintlist);
% 		
% 		hJac  = H.fhJac;
% 		
% 		Return(hJac); 
% 	);	

FreeHandJacobian(list_) 
% 	%Module({Jacmats, fhJac},
% 		
% 		% each element of list must be {DHtable_i, vars_i, gp0i } *)
% 		
% 		Jacmats = DHJacobBase(#((1)) @@ #((2)), #((3)))& /@ list;
% 		fhJac   = BlockDiag(Jacmats);
% 		
% 		Return(fhJac);
% 		
% 	);

GlobalConstraintMatrix(list_) 
% 	%Module({consmats, H},
% 		
% 		consmats = (ConstraintMatrix @@ #)& /@ list;
% 		
% 		H        = BlockDiag(consmats);
% 		
% 		Return(H);
% 		
% 	);

ConstraintMatrix(constrainttype_, axis_, Rbe_) 
% 	%Module({H, F, n,  p, RT, g, Adg},
% 		H = IdentityMatrix(6);
% 		
% 		p = {0, 0, 0};
% 		RT = Transpose(Rbe);
% 		
% 		g = RPToHomogeneous(RT, p);
% 		
% 		Adg = RigidAdjoint(g);
% 		
% 		if (constrainttype == "S", 
% 			H = Take(H, {1, 3}, {1, 6}); 
% 		);
% 		
% 		if (constrainttype == "R",  
% 		    F = Join( {0, 0, 0}, axis );
% 		    H = NullSpace({F});
% 		);
% 		
% 		if (constrainttype == "P",  
% 		    F = Join( axis, {0, 0, 0} );
% 		    H = NullSpace({F});
% 		);
% 		
% 		if (constrainttype == "C",  
% 		     H = IdentityMatrix(6);
% 		);
% 		
% 		if (constrainttype == "PC",  
% 		    H = Take(H, {3});
% 		);
% 		
% 		if (constrainttype == "PCWF",  
% 		    H = Take(H, {1, 3}, {1, 6});
% 		);
% 		
% 		if (constrainttype == "SF",  
% 		    H = H(({1, 2, 3, 6}, All));
% 		);
% 		
% 		Return(H.Adg);
% 	);

% GlobalGraspMatrix(EPlist_) 
% 	%Module({G},
% 		G = StackCols@@GraspMatrix/@EPlist;
% 		
% 		Return(G);
% 	);

GraspMatrix(EP_) 
% 	%Module({G},
% 		G = BlockMatrix(
% 			{{ IdentityMatrix(3),     ZeroMatrix(3) },
%             {           Hat(EP), IdentityMatrix(3) }}
% 		);
% 		
% 		Return(G);
% 	);
%}

%% Null stuff
%{

% only the head of the matrix has to be supplied *)
% 	
% NullBasisNumerical(A_, q0_, independent_)  
% 	%Module({m, n, Aq0},
% 		
% 		if ( MatrixQ(A(q0)) , Aq0 = A(q0), Aq0 = A@@q0);
% 		
% 		{m, n} = Dimensions( Aq0 );
% 		   eye = IdentityMatrix(n - m);
% 		   
% 		dependent = Complement(Range(n), independent);
% 		 depsel   = SelectionMatrixColumns(n, dependent);
% 		 indepsel = SelectionMatrixColumns(n, independent);
% 		   LHSmat =  Aq0.depsel;
% 		   RHSvec = -Aq0.indepsel;
% 		      res = PseudoInverse(LHSmat).RHSvec;
% 		      
% 		 For(i=1, i<=Length(independent), i++,
% 		 	
% 		 	  res = Insert(res, eye(( i )), {independent((i))} );
% 		 
% 		 );
% 		 
% 		 Return(res);
% 	);

% the whole expression A(q1,q2,...) must be supplied *)
% 	
% NullBasis(A_, independent_)  
% 	%Module({m, n, eye, dependent, depsel, indepsel, LHSmat, RHSvec, i},
% 		  
% 		   {m, n} = Dimensions(A);
% 		      eye = IdentityMatrix(n - m);
% 		dependent = Complement(Range(n), independent);
% 		 depsel   = SelectionMatrixColumns(n, dependent);
% 		 indepsel = SelectionMatrixColumns(n, independent);
% 		   LHSmat =  A.depsel;
% 		   RHSvec = -A.indepsel;
% 		      res = Inverse(LHSmat).RHSvec;
% 		      
% 		 For(i=1, i<=Length(independent), i++,
% 		 	
% 		 	  res = Insert(res, eye(( i )), {independent((i))} );
% 		 
% 		 );
% 		 
% 		 Return(res);
% 		      
% 		 
% 	);
%}

%% Rodriguez
%{
RodriguezToSpatialJac[gamma1_, gamma2_, gamma3_] :=
	Module[ {gamma, modulusgammasquared},
		gamma = {gamma1, gamma2, gamma3};
		modulusgammasquared = gamma.gamma;
		2/(1 + modulusgammasquared) { {      1,   -gamma3,    gamma2},
									  { gamma3,         1,   -gamma1},
									  {-gamma2,    gamma1,         1}
									}  	
	
	];
	

RodriguezToBodyJac(gamma1_, gamma2_, gamma3_)  
% 	%Module( {gamma, modulusgammasquared},
% 		gamma = {gamma1, gamma2, gamma3};
% 		modulusgammasquared = gamma.gamma;
% 		2/(1 + modulusgammasquared) { {      1,   gamma3,    -gamma2},
% 									  { -gamma3,         1,   gamma1},
% 									  {  gamma2,    -gamma1,         1}
% 									}  	
% 	
% 	);	

% Rodriguez parameters gamma = r tan(theta/2) *)
% RodriguezToMat(gamma1_, gamma2_, gamma3_)  
% 	%Module( {gamma, hatgamma, modulusgammasquared},
% 	
% 			gamma = {gamma1, gamma2, gamma3};
% 			hatgamma = Hat(gamma);
% 			modulusgammasquared = gamma.gamma;
% 			
% 			IdentityMatrix(3) + 2/(1 + modulusgammasquared) (hatgamma + hatgamma.hatgamma)
% 	      );

%}

%% Spline
%{


% squad(x0_,x1_,x2_,x3_,t_)  (1-2t (1 - t))*((1-t) x0 + t x3) + 2t*(1-t)((1-t)x1 + t x2)


SUBspline(p0_,p1_,p2_,p3_,t_)  
% 	Slerp(Slerp(Slerp(p0,p1,(t+2)/3.),Slerp(p1,p2,(t+1)/3.),(t+1)/2.),
%               Slerp(Slerp(p1,p2,(t+1)/3.),Slerp(p2,p3,t/3.),                       t/2.),t)
% 
% 

SBZspline(p0_,p1_,p2_,p3_,t_)  
% 	Slerp(Slerp(Slerp(p0,p1,t),Slerp(p1,p2,t),t),
%               Slerp(Slerp(p1,p2,t),Slerp(p2,p3,t),t),t)
% 

SCRspline(p0_,p1_,p2_,p3_,t_)  
% 	Slerp(Slerp(Slerp(p0,p1,t+1),Slerp(p1,p2,t),(t+1)/2.),
%               Slerp(Slerp(p1,p2,t),Slerp(p2,p3,t-1),           t/2.),         t)
% 

UBspline(p0_,p1_,p2_,p3_,t_)  
% 	Lerp(Lerp(Lerp(p0,p1,(t+2)/3.),Lerp(p1,p2,(t+1)/3.),(t+1)/2.),
%               Lerp(Lerp(p1,p2,(t+1)/3.),Lerp(p2,p3,t/3.),                       t/2.),t)
% 

BZspline(p0_,p1_,p2_,p3_,t_)  
% 	Lerp(Lerp(Lerp(p0,p1,t),Lerp(p1,p2,t),t),
%               Lerp(Lerp(p1,p2,t),Lerp(p2,p3,t),t),t)
% 
% 

CRspline(p0_,p1_,p2_,p3_,t_)  
% 	Lerp(Lerp(Lerp(p0,p1,t+1),Lerp(p1,p2,t),(t+1)/2.),
%               Lerp(Lerp(p1,p2,t),Lerp(p2,p3,t-1),           t/2.),         t)
% 

Slerp(p0_List,p1_List,t_)   %Module({costh = p0 . p1//Chop, th, sinth},
%   if (costh > 0.0,costh = Chop(costh -1.) +1.);
%   if (costh < 0.0,costh = Chop(costh +1.) -1.);
%  th = N(ArcCos(costh));
%  sinth = N(Sin(th));
% if (sinth == 0, 
%  (1-t)*p0 + t*p1,
%  (Sin(th*(1-t))/sinth)*p0 +(Sin(th*t)/sinth)*p1 ))

Lerp(p0_List,p1_List,t_)   (1-t)*p0 + t*p1

%}

%% with Spacial
%{

% QuatToBodyJac( { b0_, b1_, b2_, b3_} )  
% 	%Module({},
% 	
% 	2 * { { -b1,  b0,  b3, -b2  },
% 		  { -b2, -b3,  b0,  b1 },
% 		  { -b3,  b2, -b1,  b0 }
% 	    }
% 	
% 	);		*)
% 
% QuatToBodyJac( b0_, b1_, b2_, b3_ )  
% 	%Module({},
% 	
% 	2 * { { -b1,  b0,  b3, -b2  },
% 		  { -b2, -b3,  b0,  b1 },
% 		  { -b3,  b2, -b1,  b0 }
% 	    }
% 	); 

QuatToSpatialJac( b0_, b1_, b2_, b3_ )  
% 	%Module({},
% 	
% 	2 * { { -b1,  b0, -b3,  b2  },
% 		  { -b2,  b3,  b0, -b1 },
% 		  { -b3, -b2,  b1,  b0 }
% 	    }
% 	
% 	);
% 	

QuatToSpatialJac( { b0_, b1_, b2_, b3_} )  
% 	%Module({},
% 	
% 	2 * { { -b1,  b0, -b3,  b2  },
% 		  { -b2,  b3,  b0, -b1 },
% 		  { -b3, -b2,  b1,  b0 }
% 	    }
% 	
% 	);	*)

EulZYXToBodyJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{          -Sin(theta),         0,  1 }         
% 		 {  Cos(theta) Sin(psi),  Cos(psi),  0 }
% 		 {  Cos(theta) Cos(psi), -Sin(psi),  0 }
% 		}
% 		
% 	);

EulZYXToSpatialJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{ 0, -Sin(phi), Cos(theta) Cos(phi)  },
% 		 { 0,  Cos(phi), Cos(theta) Sin(phi)  },
% 		 { 1,        0,          -Sin(theta)  }
% 		}
% 		
% 	);

EulZYZToBodyJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{ -Cos(psi) Sin(theta), Sin(psi), 0 },
% 		 {  Sin(psi) Sin(theta), Cos(psi), 0 }
% 		 {           Cos(theta),        0, 1 }
% 		}
% 		
% 	);

EulZYZToSpatialJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{ 0, -Sin(phi), Cos(phi) Sin(theta)  },
% 		 { 0,  Cos(phi), Sin(phi) Sin(theta)  },
% 		 { 1,        0,           Cos(theta)  }
% 		}
% 		
% 	);

EulZXZToBodyJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{ Sin(theta) Sin(psi),  Cos(psi), 0 },
% 		 { Sin(theta) Cos(psi), -Sin(psi), 0 },
% 		 {          Cos(theta),         0, 1 }
% 		}
% 		
% 	);

EulZXZToSpatialJac(phi_, theta_, psi_) 
% 	%Module({},
% 		{{ 0, Cos(phi),  Sin(phi) Sin(theta) },
% 		 { 0, Sin(phi), -Cos(phi) Sin(theta) },
% 		 { 1,        0,           Cos(theta) }
% 		}
% 		
% 	);

%}

%% Dynamics 
%{

DynamicEquationsBase(DHtable_, CGtable_, Masslist_, Tensortable_, gb_, q_, qp_, v_, vp_, {Tb0_, Tne_}) 
% 	%Module({B, C, G, eqs},
% 		
% 		B = InertiaBase(DHtable, CGtable, Masslist, Tensortable, {Tb0, Tne});
% 		C = InertiaToCoriolis(B, q, qp);
% 		G = GravitationalBase(DHtable, CGtable, Masslist, gb, {Tb0, Tne});
% 		
% 		eqs = B.vp + C.v + G;
% 		
% 		Return(eqs);	
% 	);

DynamicEquations(DHtable_, CGtable_, Masslist_, Tensortable_, g0_, q_, qp_, v_, vp_) 
% 	%Module({B, C, G, eqs},
% 		
% 		B = Inertia(DHtable, CGtable, Masslist, Tensortable);
% 		C = InertiaToCoriolis(B, q, qp);
% 		G = Gravitational(DHtable, CGtable, Masslist, g0);
% 		
% 		eqs = B.vp + C.v + G;
% 		
% 		Return(eqs);	
% 	);

GravitationalBase(DHtable_, CGtable_, Masslist_, g0_, {Tb0_, Tne_}) 
% 	%Module({Gvec, k, n = Length(DHtable)},
% 		
% 		Gvec = (ZeroMatrix(1, n))((1));
% 		
% 		For( k = 1, k <= n, k++,
% 			
% 			R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
% 		    Jk  =  CGJacobBaseDyn(DHtable, CGtable, {Tb0, Tne}, k);
% 		    
% 		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
% 		    
% 		    Jvk =  Jk((1;;3, All));
% 		    % Jok =  Jk((4;;6, All)); *)
% 		    
% 		    JvkT = Transpose(Jvk);
% 		    % JokT = Transpose(Jok);  *)
% 		    
% 		    Gvec += -Masslist((k)) JvkT.g0; 
% 		);
% 		Return(Gvec);
% 		
% 	);

Gravitational(DHtable_, CGtable_, Masslist_, g0_) 
% 	%Module({Gvec, k, n = Length(DHtable)},
% 		
% 		Gvec = (ZeroMatrix(1, n))((1));
% 		
% 		For( k = 1, k <= n, k++,
% 			
% 			R0k  = RigidOrientation( DHFKine(DHtable, k) );
% 		    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
% 		    
% 		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
% 		    
% 		    Jvk =  Jk((1;;3, All));
% 		    % Jok =  Jk((4;;6, All)); *)
% 		    
% 		    JvkT = Transpose(Jvk);
% 		    % JokT = Transpose(Jok);  *)
% 		    
% 		    Gvec += -Masslist((k)) JvkT.g0; 
% 		);
% 		Return(Gvec);
% 		
% 	);

InertiaToCoriolisNotChristoffel(M_, q_, qp_)  
% 	%Module(
% 		{n = Length(M), Cmat, i, j, k},
% 		Cmat = ZeroMatrix(n);
% 		
% 	   	   
% 	   	   For(i = 1, i <= n, ++i,
%              For(j = 1, j <= n, ++j,
%                 Cmat((i,j)) = Sum(
%            	                  ( D( M((i, j)), q((k)) ) - (1/2) D( M((j, k)), q((i)) ) ) qp((k)), {k, 1, n}
%            	                     )
%              )
% 	   	   );
%            Return(Cmat);	                  
% 	);

InertiaToCoriolis(M_, q_, qp_)  
%   %Module(
%     {Cmat, i, j, k, n = Length(M)},
% 
%     % Brute force calculation *)
%     Cmat = ZeroMatrix(n);
% 
%     For(i = 1, i <= n, ++i,
%       For(j = 1, j <= n, ++j,
%         For(k = 1, k <= n, ++k,
%           Cmat((i,j)) += (1/2) * qp((k)) * (D(M((i,j)), q((k))) + D(M((i,k)), q((j))) - D(M((j,k)), q((i))))
% 	       )
%          )
%        );
%    Cmat
%   );

InertiaBase(DHtable_, CGtable_, Masslist_, Tensortable_, {Tb0_, Tne_}) 
% 	%Module(
% 		{Bmat,
% 			k, n = Length(DHtable),
% 		 R0k, Jk, Jvk, Jok, JvkT, JokT},
% 		 
% 		 Bmat = ZeroMatrix(n);
% 		
% 		For( k = 1, k <= n, k++,
% 			
% 			R0k  = RigidOrientation( DHFKine(DHtable, {Tb0, Tne}, k) );
% 		    Jk  =  CGJacobBaseDyn(DHtable, CGtable, {Tb0, Tne}, k);
% 		    
% 		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
% 		    
% 		    Jvk =  Jk((1;;3, All));
% 		    Jok =  Jk((4;;6, All));
% 		    
% 		    JvkT = Transpose(Jvk);
% 		    JokT = Transpose(Jok);
% 		    
% 		    Bmat += Masslist((k)) JvkT.Jvk + JokT.R0k.Tensortable((k)).Transpose(R0k).Jok; 
% 		);
% 		Return(Bmat);
% 	);

Inertia(DHtable_, CGtable_, Masslist_, Tensortable_) 
% 	%Module(
% 		{Bmat,
% 			k, n = Length(DHtable),
% 		 R0k, Jk, Jvk, Jok, JvkT, JokT},
% 		 
% 		 Bmat = ZeroMatrix(n);
% 		
% 		For( k = 1, k <= n, k++,
% 			
% 			R0k  = RigidOrientation( DHFKine(DHtable, k) );
% 		    Jk  =  CGJacob0Dyn(DHtable, CGtable, k);
% 		    
% 		    % attenzione: questo Jacobiano \(EGrave) relativo al baricentro *)
% 		    
% 		    Jvk =  Jk((1;;3, All));
% 		    Jok =  Jk((4;;6, All));
% 		    
% 		    JvkT = Transpose(Jvk);
% 		    JokT = Transpose(Jok);
% 		    
% 		    Bmat += Masslist((k)) JvkT.Jvk + JokT.R0k.Tensortable((k)).Transpose(R0k).Jok; 
% 		);
% 		Return(Bmat);
% 	);

Regressor(DHtable_, q_, qp_, v_, vp_, t_, g0_) 
% 	%Module({n, Yk, Y},
% 		
% 		n = Length(DHtable);
% 		Y = Regressor(DHtable, q, qp, v, vp, t, g0, 1);
% 		For( k = 2, k <= n, k++,
% 			
% 			Yk = Regressor(DHtable, q, qp, v, vp, t, g0, k);
% 			Y  = StackCols(Y, Yk);
% 		);
% 		
% 		Return(Y);
% 	);

Regressor[DHtable_, q_, qp_, v_, vp_, t_, g0_, k_]:=
	Module[{X0p, X1a, X1b, X1ap, X1bp, X1p, X2p, par,
		    W0, W1, W2,
		    Z0, Z1,
		    R0k, R0kT,
		    Jk, Jvk, Jok, JvkT, JokT,
		    Jvkp, JvkpT, Jokp, JokpT,
		    JvTJv,
		    T1, T2, T3, TT, 
		    E1, E2, E3, E4, E5, E6, EE,
		    Y0, Y1, Y2, Y},
		    
		    R0k  = RigidOrientation[ DHFKine[DHtable, k] ];
		    R0kT = Transpose[R0k];
		    Jk   =  DHJacob0Dyn[DHtable, k];
		    Jvk  =  Jk[[1;;3, All]];
		    Jok  =  Jk[[4;;6, All]];
		    
		    JvkT = Transpose[Jvk];
		    JokT = Transpose[Jok];
		    JvTJv = JvkT.Jvk;
		    
		    Jvkp  = TensorDerivative[Jvk, q].qp;
		    JvkpT = Transpose[Jvkp];
		    
		    Jokp  = TensorDerivative[Jok, q].qp;
		    JokpT = Transpose[Jokp]; 
		    
		    (* termini da d   dT1                                        *)
		    (*            -   -                                          *)
		    (*            dt  dqp         con T1 = v.B(q).qp             *)
		    
		    
 			X0p = Transpose[{(TensorDerivative[JvTJv, q].qp).v  + JvTJv.vp}];
 			
 			
 			(* X1p = (   JvkpT.Hat[Jok.v] + JvkT.(Hat[Jokp.v] + Hat[Jok.vp])
 				   - JokpT.Hat[Jvk.v]  - JokT.(Hat[Jvkp.v] + Hat[Jvk.vp]) ).R0k +
 				      (JvkT.Hat[Jok.v] - JokT.Hat[Jvk.v]).(D[R0k, t]); *)
 			
 			T1 = SparseArray[{{2, 3} -> -1,  {3, 2} ->  1, {3,3} -> 0}];
 			T2 = SparseArray[{{1, 3} ->  1,  {3, 1} -> -1, {3,3} -> 0}];
 			T3 = SparseArray[{{1, 2} -> -1,  {2, 1} ->  1, {3,3} -> 0}];
 			
 			TT = {T1, T2, T3};
 			TT = Transpose[TT, {2,1,3}];
 			
 			X1a  = JokT.R0k.TT.R0kT.Jvk;
 			
 			X1ap = Transpose[Table[(TensorDerivative[X1a[[All,i,All]], q].qp).v, {i,1,3}]];
 			
 			X1b  = -JvkT.R0k.TT.R0kT.Jok;  
 			
 			X1bp = Transpose[Table[(TensorDerivative[X1b[[All,i,All]], q].qp).v, {i,1,3}]];
 			
 			
 			X1p = (X1a + X1b).vp + X1ap + X1bp;
 			
 			E1 = SparseArray[{{1, 1} -> 1,  {3, 3} -> 0}];
 			E2 = SparseArray[{{1, 2} -> -1, {2, 1} -> -1 , {3, 3} -> 0}];
 			E3 = SparseArray[{{1, 3} -> -1, {3, 1} -> -1 , {3, 3} -> 0}];
 			E4 = SparseArray[{{2, 2} -> 1,  {3, 3} -> 0}];
 			E5 = SparseArray[{{2, 3} -> -1, {3, 2} -> -1 , {3, 3} -> 0}];
 			E6 = SparseArray[{{3, 3} -> 1}];
 			EE = {E1, E2, E3, E4, E5, E6};
 			EE = Transpose[EE, {2,1,3}];
 			
 			par = JokT.R0k.EE.R0kT.Jok;
 			
 			X2p = Transpose[Table[(TensorDerivative[par[[All,i,All]], q].qp).v, {i,1,6}]] + par.vp;
 		
 			
 			(* termini da dT2                                    *)
		    (*            -                                      *)
		    (*            dq           con T2 = (1/2)v.B(q).qp   *)
		    
		    W0 = (1/2) D[ v.JvkT.Jvk.qp ,{q}];
		    
		    W1 = (1/2) Transpose[ D[ (v.JvkT.Hat[Jok.qp] - v.JokT.Hat[Jvk.qp]).R0k  ,{q}] ];
 			
 			W2 = (1/2) Transpose[ D[ v.JokT.R0k.EE.Transpose[R0k].Jok.qp ,{q}] ];
 			
 			
 			(* termini da d U   *)
 			(*            -     *)
 			(*            d q   *)
 			
 			Z0 = -JvkT.g0;
 			
 			(*
 			If[ k == 1,
 			Z1 = -Transpose[ { D[ g0.R0k , {q}] } ],
 			Z1 = -Transpose[   D[ g0.R0k , {q}]   ]
 			 ];
 			 *)
 			 
 			Z1 = -Transpose[   D[ g0.R0k , {q}]   ]; 
 			
 			Y0 = X0p - W0  + Z0;
 			Y1 = X1p - W1  + Z1;
 			Y2 = X2p - W2;
 			
 			Y  = StackCols[Y0, Y1, Y2];
 			Return[Y]; 
 			 
 			 
 			 (* test:   Return[X2p]; *)	
	]; 

TensorDerivative(J_, q_) 	
% 	%Module({DJ, res },
% 		DJ  = D(J, {q});
% 		res =(1/2) (DJ + Transpose(DJ, {1, 3, 2}));
% 		Return(res);
% 	);

%}

%%
%{
Mio                                         Grabovic

        (RigidAdjoint  InverseRigidAdjoint)
NA                                          adjoint
NA                                          ad
NA                                          AD_FWKin
NA                                          adjointInv
NA                                          adjointStar
NA                                          adStar

NA                                          axisAngleToQuat
        (BodyJacobian )
NA                                          BodyJac
NA                                          BodyJacAD
NA                                          BodyJacDerivative

[round]                                     chop

NA                                          diffRotZ2D

SkewExp                                     expSkew
TwistExp                                    expTw
        (ForwardKinematics)
NA                                          FWKin

Hat                                     !=  hat
Skew                                        NA
AxisToSkew                                  NA
HatInv                                      vec
Unskew                                      NA
SkewToAxis                                  NA

NA                                          hom_mat
RigidInverse                                hom_mat_inv
NA                                          hom_mat_inv2D

NA                                          omegaFromQuat
NA                                          omegaFromQuat2

QuatToMat                                          quatToRot

RigidInverse                                rigidInverse

NA                                          RK4_integrator

NA                                          rotNtheta
NA                                          rotNthetaQuat

NA                                          rotToAxisAngle

RotX                                        rotX
NA                                          rotXvectorized
RotY                                        rotY
RotZ                                        rotZ
NA                                          rotZ2D

NA                                          spatialJac
NA                                          spatialJacDerivative

NA                                          spatiaTwistDerivative

HomX                                        TrotX
HomY                                        TrotY
HomZ                                        TrotZ
HomX                                        Ttx
HomY                                        Tty
HomZ                                        Ttz

NA                                          unitTwist