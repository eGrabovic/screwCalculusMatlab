function mat = QuatToMat(quat)
%{
 QuatToMat(q_List)  
Module({q0= q((1)), q1=q((2)), q2 = q((3)), q3 = q((4))},0
   %Module({d23 = 2 q2 q3, d1 = 2 q0 q1,
          d31 = 2 q3 q1, d2 = 2 q0 q2,
          d12 = 2 q1 q2, d3 = 2 q0 q3,
          q0sq = q0^2, q1sq = q1^2, q2sq = q2^2, q3sq = q3^2},
          {{q0sq + q1sq - q2sq - q3sq, d12 - d3, d31 + d2},
           {d12 + d3, q0sq - q1sq + q2sq - q3sq, d23 - d1},
           {d31 - d2, d23 + d1, q0sq - q1sq - q2sq + q3sq}} )); 
%}
% TODO: check
    scalar = quat(1);
    vectorial = quat(2:end);
    
    Im = eye(3);
    vect_hat = Skew(vectorial);

    mat = Im + 2*vect_hat* (scalar * Im + vect_hat);
end