function [u, debug] = u_NMPC(X)
x0=X(1); % x position (obsolete)
y0=X(2); % y position (obsolete)
v0=(X(3)); % velocity
beta0=X(4); % side slip angle
psi0=X(5); % yaw angle
omega0=X(6); % yaw rate
%x_dot=X(7); % longitudinal velocity (obsolete)
%y_dot=X(8); % lateral velocity (obsolete)
psi_dot=X(9); % yaw rate (redundant)
varphi_dot=(X(10)); % wheel rotary frequency

% load("constraints", "ref_func", "x0_f", "lb_f", "ub_f");
load("constraints", "k_ref_func", "x0_f", "lb", "ub");

p0 = [x0; y0];

s = find_s(p0);
[k_ref, psi0_ref, x0_ref, y0_ref] = referencePath(s);
e_y_0 = cos(psi0_ref) * (y0 - y0_ref) - sin(psi0_ref) * (x0 - x0_ref);
e_psi_0 = psi0 - psi0_ref;
t0 = 0;
q0 = [e_y_0; e_psi_0; v0; beta0; omega0; t0];
x0 = x0_f(q0);
ref = k_ref_func(s);
% lb = lb_f(q0);
% ub = ub_f(q0);
f = @costFunction;
nonlin = @(x) c_nonlin(x, q0, ref);
opts = optimoptions('fmincon','Algorithm','sqp', 'Display','none', ...
    "MaxFunctionEvaluations", 3000, "UseParallel",true);
sol = fmincon(f, x0, [], [], [], [], lb, ub, nonlin, opts);
u1 = sol(1);
delta = sol(2);

load("T_inv.mat", "T_inv");
if u1 >= 0
    R = 0.302;
    i0 = 3.91;
    i_g=[3.91 2.002 1.33 1 0.805];
    % T_want = (1/R * i0 * i_g(G))^-1 * u(1);
    [u1_max, G] = u1_max_f(v0);
    u1_real = min(u1, u1_max);
    T_real = (1/R * i0 * i_g(G))^-1 * u1_real;
    n = v0 * i_g(G) * i0 / R;
    phi = T_inv(T_real, n);
    Fb = 0;
else
    G = 1; % Does not matter
    phi = 0;
    Fb = -u1;
end

zeta = 1;
% u = [phi, G, Fb, zeta, delta];
u = [delta, G, Fb, zeta, phi];
debug = [e_y_0, e_psi_0, s, k_ref];