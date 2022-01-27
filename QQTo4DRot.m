function rot = QQTo4DRot(p, q)
% QQTO4DROT(p, q) returns the 4D rotation matrix ROT from two unit 
%   quaternions P and Q

  q0 = q(1);
  q1 = q(2);
  q2 = q(3);
  q3 = q(4);
  p0 = p(1);
  p1 = p(2);
  p2 = p(3);
  p3 = p(4);


rot = [
        [p0*q0 + p1*q1 + p2*q2 + p3*q3, ...
         p1*q0 - p0*q1 - p3*q2 + p2*q3, ...
         p2*q0 + p3*q1 - p0*q2 - p1*q3, ...
         p3*q0 - p2*q1 + p1*q2 - p0*q3
        ];
        [-(p1*q0) + p0*q1 - p3*q2 + p2*q3, ... 
         p0*q0 + p1*q1 - p2*q2 - p3*q3, ...
         -(p3*q0) + p2*q1 + p1*q2 - p0*q3, ...
         p2*q0 + p3*q1 + p0*q2 + p1*q3
        ];
        [-(p2*q0) + p3*q1 + p0*q2 - p1*q3, ...
         p3*q0 + p2*q1 + p1*q2 + p0*q3, ...
         p0*q0 - p1*q1 + p2*q2 - p3*q3, ...
         -(p1*q0) - p0*q1 + p3*q2 + p2*q3
        ];
        [-(p3*q0) - p2*q1 + p1*q2 + p0*q3, ...
         -(p2*q0) + p3*q1 - p0*q2 + p1*q3, ...
         p1*q0 + p0*q1 + p3*q2 + p2*q3, ...
         p0*q0 - p1*q1 - p2*q2 + p3*q3
        ]
      ];
end