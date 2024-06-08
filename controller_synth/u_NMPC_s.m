function u = u_NMPC_s(s, q0)
addpath("..\pathGen\")
% load("constraints", "ref_func", "x0_f", "lb_f", "ub_f");
load("constraints", "k_ref_func", "x0_f", "lb", "ub");
k_ref = k_ref_func(s);
q0(6) = 0;
x0 = x0_f(q0);
% lb = lb_f(q0);
% ub = ub_f(q0);
f = @costFunction;
nonlin = @(x) c_nonlin(x, q0, k_ref);
opts = optimoptions('fmincon','Algorithm','sqp', 'Display','none', ...
    "MaxFunctionEvaluations", 3000, "UseParallel",true);
sol = fmincon(f, x0, [], [], [], [], lb, ub, nonlin, opts);
% u = zeros(4, 1);
% u([1,2,4]) = sol(1:3);
u = sol(1:2);
