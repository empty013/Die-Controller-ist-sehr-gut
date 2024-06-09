% function synth_NMPC(s, q0, ref_f)
clear;
% q0 = [0; 0; 5; 0; 0; 0];


addpath("..\pathGen\")
ref_f = @referencePath;
load("f_s_sym.mat", "f_sym");
syms e_y e_psi v beta omega t
% syms u1 F_b zeta delta
% syms u1 F_b delta
syms u1 delta
syms k_ref psi_ref
q = [e_y; e_psi; v; beta; omega; t];
% u = [u1; F_b; zeta; delta];
% u = [u1; F_b; delta];
u = [u1; delta];
ref = [k_ref];
syms q0 [length(q), 1]
syms s
% s = 0;


% Setup prediction horizon
v_ref = 10;
s_end = 0.4;
s_delta = 0.02;
N = s_end / s_delta;
s_arr = linspace(s, s + s_end, N);
% syms q_arr [length(q), N]
syms q_arr [length(q), N-1]
syms u_arr [length(u), N-1]
syms k_ref_arr [1, N]
x = [];
x0 = [];
A = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
c_eq = [];
c_ineq = [];
Q_tot = [];
R_tot = [];
q_ref = [];

% q_i = q_arr(:, 1);
q_i = q0;
q_ref_i = [0; 0; v_ref; 0; 0; 0];
% ref_i = @(s) referencePath(s);
% x = [x; q_i];
% x0 = [x0; q0];
% lb = [lb; q0];
% ub = [ub; q0];
% Q_tot = blkdiag(Q_tot, Q);
% t_cost = 2*s_end;
Q = diag([5e1; 1e-3; 1e-4; 1e-1; 1e1; 1e-2]);
e_psi_beta_mixed_weight = 1e0;
Q(2, 2) = Q(2,2) + e_psi_beta_mixed_weight;
Q(4, 4) = Q(4,4) + e_psi_beta_mixed_weight;
Q(2, 4) = Q(2,4) - e_psi_beta_mixed_weight;
Q(4, 2) = Q(4,2) - e_psi_beta_mixed_weight;
R = diag([0; 1e-1]);
for ii = 1:N-1
    if mod(ii, 20) == 0
        ii
    end
    u_i = u_arr(:, ii);
    x = [x; u_i];
    % x0 = [x0; [5000; 0; 0.5; 0]];
    % lb = [lb; [0; 0; 0; -0.53]];
    % ub = [ub; [10000; 15000; 1; 0.53]];
    x0 = [x0; [0; 0]];
    lb = [lb; [-1.5; -0.53]];
    ub = [ub; [1; 0.53]];
    R_tot = blkdiag(R_tot, R);
    
    
    % ref_i(1:2,1) = ref_f(s_arr(ii));
    k_ref_i = k_ref_arr(ii);
    f_i = subs(f_sym, [q; u; ref], [q_i; u_i; k_ref_i]);
    q_i_end = q_i + s_delta * f_i;
    % k1 = subs(f_sym, [q; u; ref], [q_i; u_i; ref_i]);
    % k2 = subs(f_sym, [q; u; ref], [q_i + s_delta/2 * k1; u_i; ref_i]);
    % k3 = subs(f_sym, [q; u; ref], [q_i + s_delta/2 * k2; u_i; ref_i]);
    % k4 = subs(f_sym, [q; u; ref], [q_i + s_delta * k3; u_i; ref_i]);
    % q_i_end = simplify(q_i + s_delta/6 * (k1 + 2*k2 + 2*k3 + k4));
    
    % q_i = q_arr(:, ii+1);
    q_i = q_arr(:, ii);

    q_ref = [q_ref; q_ref_i];
    x = [x; q_i];
    x0 = [x0; q0];
    lb = [lb; [-Inf; -Inf; -Inf; -Inf; -Inf; -Inf]];
    ub = [ub; [Inf; Inf; Inf; Inf; Inf; Inf]];
    % lb = [lb; [-0.2; -1; 0; -Inf; -Inf; -Inf]];
    % ub = [ub; [0.2; 1; Inf; Inf; Inf; Inf]];
    % lb = [lb; [0; -1; 0; -Inf; -Inf; -Inf]];
    % ub = [ub; [0; 1; Inf; Inf; Inf; Inf]];

    c_eq = [c_eq; [q_i_end - q_i]];
    Q_tot = blkdiag(Q_tot, Q);
end

% lb(end-5) = 0;
% ub(end-5) = 0;
% lb(end-5) = 0;
% ub(end-5) = 0;
J = (q_arr(:) - q_ref).' * Q_tot * (q_arr(:) - q_ref) ...
    + u_arr(:).' * R_tot * u_arr(:) ...
    + 1e3 * q_i(end) + 0e4 * (q_i(2) - q_i(4))^2 + 1e5 * q_i(1)^2;

s_arr_f = matlabFunction(s_arr, "Vars", s);
k_ref_func = @(s) referencePathCurv(s_arr_f(s));
J_f = matlabFunction(J, "Vars", {x}, "File", "costFunction");
x0_f = matlabFunction(x0, "Vars", {q0});
% lb_f = matlabFunction(lb, "Vars", {q0});
% ub_f = matlabFunction(ub, "Vars", {q0});
tic
c_eq_f = matlabFunction(c_eq, "Vars", {x; q0; k_ref_arr}, "File", "c_eq_f");
toc


% c_eq_dx = jacobian(c_eq, x);
% e_eq_dx_sparse = sparse(c_eq_dx);
% c_eq_dx_f = matlabFunction(c_eq_dx, ...
%       "Vars", {x; q0; ref_arr}, "File", "c_eq_dx_f");

% save("constraints", "ref_func", "x0_f", "lb_f", "ub_f");
save("constraints", "k_ref_func", "x0_f", "lb", "ub", "s_end", "s_delta");
