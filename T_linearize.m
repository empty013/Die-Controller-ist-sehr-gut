clear; close all
syms v phi g n
% R = 0.302;
% i0 = 3.91;
% i_g=[3.91 2.002 1.33 1 0.805];
% G = 2;

% Higher values of v(=19) for 2nd gear not recommended, check plot
n_max = 480;
n_res = 100;
phi_res = 100;
n_arr = linspace(0, n_max, n_res);
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

%%
% Check quality of phi = T_inv(T_M)
figure(2) 
hold on
view(3)
xlabel("\phi_{true}")
ylabel("n")
zlabel("\phi_{inv}")
title("Check linearization")

figure(3)
hold on
view(3)
xlabel('T_M')
ylabel('n')
zlabel('\phi')
title("\phi = T^{-1}(T, n)")

% figure(3)
% xlabel("\phi")
% ylabel("v")
% zlabel("T_M")
% hold on

options_fsolve = optimoptions(@fsolve, "Display", "iter", ...
    "Algorithm", "trust-region-dogleg");


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
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'off' );

coeffs2 = num2cell(coeffvalues(fitresult));
coeffs = {0.232540660594679	0.114393289007004	0.0779463002016656	...
    0.0115624070879991	-0.00427974336827419	-0.0102118361339179	...
    -0.00145644719122835	0.0184573885301343	0.0300930337094600	...
    0.00866855693966725	-0.00650393473156937	0.0309654815549738	...
    0.0440918543444869	0.0754763843366104	0.0494180156914597	...
    0.00861505668499663	0.00347726441814111	0.0148511673968274	...
    0.0162403456177435	0.0308870560261451	0.0231914355019211};
T_inv = sfit(ft, coeffs2{:});

figure(2)
n_idx = 2;
% surf(T_sol_mesh, n_sol_mesh, fitresult(T_sol_mesh, n_sol_mesh));
% phi_arr_test = phi_arr(1:ind(n_idx));
% n_arr_test = n_mesh(n_idx, 1:ind(n_idx));
% T_arr_test = T_mesh(n_idx, 1:ind(n_idx));
% phi_inv_test = T_inv(T_arr_test, n_arr_test);
% plot3(phi_arr_test, n_arr_test, phi_inv_test);
% plot3(phi_arr_test([1,end]), n_arr_test([1,end]), ...
%     phi_arr_test([1,end]), 'r')

plot3(phi_sol_mesh.', n_sol_mesh.', T_inv(T_sol_mesh, n_sol_mesh).')
plot3(phi_sol_mesh(:,[1,end]).', n_sol_mesh(:,[1,end]).', ...
    phi_sol_mesh(:,[1,end]).', 'r')