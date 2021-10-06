syms beta  a b phi gamma Ix Iy m phi_dot gamma_dot beta_dot beta_dotdot
assume(beta, 'real')
assume(gamma, 'real')
assume(phi, 'real')


beta = ADvar(ADvar(beta,beta_dot), ADvar(beta_dot, beta_dotdot));
phi = ADvar(ADvar(phi,phi_dot), ADvar(phi_dot, 0));
gamma = ADvar(ADvar(gamma, gamma_dot), ADvar(gamma_dot, 0));

Rgq = rotY(beta);
dgq = [0; a+b;0];
Iq = eye(3).*[Ix;Iy;Iy];
Ig = Rgq*Iq*Rgq.' + m.*sC.hat(dgq)*sC.hat(dgq).';