clear; close all
syms v phi g n
R = 0.302;
i0 = 3.91;
i_g=[3.91 2.002 1.33 1 0.805];
G = 2;

% Higher values of v(=19) for 2nd gear not recommended, check plot
n_res = 100;
phi_res = 100;
n_arr = linspace(0, 480, n_res);
phi_arr = linspace(0, 1, phi_res);
[phi_mesh, n_mesh]  = meshgrid(phi_arr, n_arr);


% n = v * i_g(G) * i0 / R;
T_M=200*phi*(15-14*phi)- ...
    200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); 
T_M_f = matlabFunction(T_M, "Vars", [phi, n]);

figure(1);
xlabel("\phi")
ylabel("n")
zlabel("T_M")
hold on

h = surf(phi_mesh, n_mesh, T_M_f(phi_mesh, n_mesh));
set(h,'LineStyle','none')
T_mesh = T_M_f(phi_mesh, n_mesh);
[T_max, ind] = max(T_mesh, [], 2); % maximum w.r.t. v
plot3(phi_arr(ind), n_arr, T_max, 'r-', 'LineWidth',2)

% for g = 1:1
%     % n = v * i_g(g) * i0 / R;
%     % T_M=200*phi*(15-14*phi)- ...
%         % 200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); % motor torque
%     % fplot(T_M, [0, 5])
%     surf(phi_mesh, v_mesh, T_M_f(phi_mesh, v_mesh));
%     tmp1 = phi_arr(ind);
%     tmp2 = v_arr;
%     plot3(tmp1, tmp2, m_, 'r-', 'LineWidth',2)
% end

%%
% Check quality of phi = T_inv(T_M)
figure(2) 
hold on
view(3)
xlabel("\phi_{true}")
ylabel("n")
zlabel("\phi_{inv}")

figure(3)
hold on
view(3)
xlabel('T_M')
ylabel('n')
zlabel('\phi')

% figure(3)
% xlabel("\phi")
% ylabel("v")
% zlabel("T_M")
% hold on

options_LSQ = optimoptions(@lsqnonlin, "Display","none");
options_fsolve = optimoptions(@fsolve, "Display", "iter", ...
    "Algorithm", "trust-region-dogleg");
func_model = @T_inv_poly;
param_num = 7;
syms c [1 param_num]

% n_max = 500
T_n = 50;
n_n = 50;
n_step = floor(n_res / n_n);
n_indices = 1:n_step:n_res;

phi_sol_mesh = zeros(n_n, T_n);
% v_arr_true = linspace(0, 19, v_n);
T_sol_mesh = zeros(n_n, T_n);
n_sol_mesh = zeros(n_n, T_n);
% [v_mesh_true, T_mesh_true] = meshgrid()
for id = 1:length(n_indices)
    n_idx = n_indices(id);
    phi_max_idx = ind(n_idx);
    phi_max = phi_arr(phi_max_idx);
    
    % phi_arr_true = linspace(0, phi_max, T_n);
    n_arr_true = repmat(n_mesh(n_idx, 1), 1, T_n);
    T_min = T_mesh(n_idx, 1);
    T_max = T_mesh(n_idx, phi_max_idx);
    T_arr_true = linspace(T_min, T_max, T_n);
    % T_arr_true = T_M_f(phi_arr_true, v_arr_true);

    T_res = @(phi) T_M_f(phi, n_arr_true) - T_arr_true;
    % phi0 = zeros(size(T_arr_true)) + 0.2;
    phi0 = linspace(0, phi_max, T_n);
    [phi_sol, fval] = fsolve(T_res, phi0, options_fsolve);
    T_sol_mesh(id, :) = T_arr_true;
    n_sol_mesh(id, :) = n_arr_true;
    phi_sol_mesh(id, :) = phi_sol;

    % T_inv_res_func = @(params) res_f(params, T_arr_true, phi_sol, ...
    %                                     func_model);
    % T_res_func = @(params) res_f(params, phi_arr_true, T_arr_true);
    % params0 = [1,1,1,1,1];
    % params0 = zeros(1, param_num) +0.1;
    % params0 = randn(1, param_num);
    % [params_sol, res_norm] = lsqnonlin(T_inv_res_func, params0, ...
    %                             [], [], options_LSQ);
    % T_inv = @(T) func_model(T, params_sol);

    figure(3)
    plot3(T_arr_true, n_arr_true, phi_sol);
    % plot3(T_arr_true, v_arr_true, T_inv(T_arr_true))
    
    % figure(1)
    % plot3(phi_arr_true, v_arr_true, ...
    %     T_approx(phi_arr_true, params_sol), 'g', 'LineWidth', 2)
    
    % figure(2)
    % T_inv_subs = simplify(subs(T_inv_sym, c, params_sol));
    % T_inv = matlabFunction(T_inv_subs);
    % phi_arr_test = phi_arr(1:ind(v_idx));
    % v_arr_test = v_mesh(v_idx, 1:ind(v_idx));
    % T_arr_test = T_mesh(v_idx, 1:ind(v_idx));
    % phi_inv_test = T_inv(T_arr_test);
    % plot3(phi_arr_test, v_arr_test, phi_inv_test);
    % plot3(phi_arr_test([1,end]), v_arr_test([1,end]), ...
    %     phi_arr_test([1,end]), 'r')
    
    % eval_abs = abs(phis_inv - phis(1:ind(idx)));
    % eval_data = eval_abs ./ phis(1:ind(idx));
    % [~, eval_ind] = max(eval_data(2:end));

    % tmp = T_approx(phis_inv, params_sol)
    % figure(3)
    % plot3(phi_arr(1:ind(idx)), v_mesh(idx, 1:ind(idx)), ...
    %     T_approx(phis_inv, params_sol), 'k', 'LineWidth', 2)

    % figure(1)
    % plot3(phi_arr_true, v_arr_true, ...
    %     T_inv(T_mesh(idx, 1:ind(idx))));
    % plot3(phi_arr(1:ind(idx)), v_arr_true, ...
    %     T_approx(phis_inv, params_sol), 'k', 'LineWidth', 2)
end

[xData, yData, zData] = prepareSurfaceData( T_sol_mesh, n_sol_mesh, phi_sol_mesh );

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );



function res = res_f(params, x_target, y_target, func)
    % res = f(x) - y != 0
    res = func(x_target, params) - y_target;
end

function func = T_inv_poly(x, params)
    x = x/800;
    func = ... % params(1) ...
        + params(1) .* x .* exp(params(2).*(x/1)) ...
        + params(3) .* x.^2 ...
        + params(4) .* x ...
        + params(5) .* exp(params(6) .* (x/1) + params(7));
        % + params(4) .* x.^3; 
        % + params(4) .* x.^4 ...
        % + params(4) .* x.^6;
end

function func = T_inv_sin(x, params)
    func = params(1) ...
        + params(2) .* cos(params(3) .* x) ...
        + params(4) .* sin(params(3) .* x) ...
        + params(5) .* cos(2 * params(3) .* x) ...
        + params(6) .* sin(2 * params(3) .* x);
end

function func = T_inv_approx(x, params)
    func = params(1) ... 
        + params(2) .* sinh(params(3) .* x) ...
        + params(4) .* x ...
        + params(5) .* x.^2 ...
        + params(6) .* cosh(params(7) .* x);
end