function out = RK4_integrator(x, u, xdot, dt, Nsteps)
%
%
%

h = dt./Nsteps;

x_end = x;

for i = 1:Nsteps
    
    v_1 = xdot(x_end, u);
    v_2 = xdot(x_end + 0.5.*h.*v_1, u);
    v_3 = xdot(x_end + 0.5.*h.*v_2, u);
    v_4 = xdot(x_end + v_3.*h, u);
    x_end = x_end + (1/6).*(v_1 + 2.*v_2 + 2.*v_3 + v_4).*h;

end
    
out = x_end;
end