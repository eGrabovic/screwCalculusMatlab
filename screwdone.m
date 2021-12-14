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

%}  