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
[T_max_arr, ind] = max(T_mesh, [], 2); % maximum w.r.t. v
plot3(phi_arr(ind), n_arr, T_max_arr, 'r-', 'LineWidth',2)

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

options_fsolve = optimoptions(@fsolve, "Display", "iter", ...
    "Algorithm", "trust-region-dogleg");
T_n = 10;
n_n = 10;
n_step = floor(n_res / n_n);
n_indices = 1:n_step:n_res;

phi_sol_mesh = zeros(n_n, T_n);
phi_pred_mesh = zeros(n_n, T_n);
T_pred_mesh = zeros(n_n, T_n);
n_pred_mesh = zeros(n_n, T_n);
res_sum = 0; % Check for quality of T^-1, should be small
% Numerically calculate phi = T^-1(T, n) and collect in mesh for fitting.
% Done row-wise for each n in. 
for id = 1:length(n_indices)
    n_idx = n_indices(id);
    n_arr_true = repmat(n_mesh(n_idx, 1), 1, T_n);
    phi_max_idx = ind(n_idx);
    phi_max = phi_arr(phi_max_idx);
    
    T_min = T_mesh(n_idx, 1);
    T_max = T_mesh(n_idx, phi_max_idx);
    T_arr_true = linspace(T_min, T_max, T_n);

    T_res = @(phi) T_M_f(phi, n_arr_true) - T_arr_true;
    phi0 = linspace(0, phi_max, T_n);
    phi_pred_mesh(id, :) = phi0;
    [phi_sol, fval] = fsolve(T_res, phi0, options_fsolve);
    res_sum = res_sum + sum(fval);
    T_pred_mesh(id, :) = T_arr_true;
    n_pred_mesh(id, :) = n_arr_true;
    phi_sol_mesh(id, :) = phi_sol;

    figure(3)
    plot3(T_arr_true, n_arr_true, phi_sol);
end

[xData, yData, zData] = prepareSurfaceData( T_pred_mesh, n_pred_mesh, ...
                                            phi_sol_mesh );

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, ~, output] = fit( [xData, yData], zData, ft, 'Normalize', 'off' );

coeffs2 = num2cell(coeffvalues(fitresult));
% For later, this way can create 2d poly from vector.
% ATTENTION: Curve fitting toolbox needed!!!!!!!!!!!!!!!!??????????????
coeffs = {0.232540660594679	0.114393289007004	0.0779463002016656	...
    0.0115624070879991	-0.00427974336827419	-0.0102118361339179	...
    -0.00145644719122835	0.0184573885301343	0.0300930337094600	...
    0.00866855693966725	-0.00650393473156937	0.0309654815549738	...
    0.0440918543444869	0.0754763843366104	0.0494180156914597	...
    0.00861505668499663	0.00347726441814111	0.0148511673968274	...
    0.0162403456177435	0.0308870560261451	0.0231914355019211};
T_inv = sfit(ft, coeffs2{:});
save("..\template_Matlab\T_inv", "T_inv");

figure(2)
% For plotting only the selected values of the mesh that were also fitted.
% plot3(phi_sol_mesh.', n_sol_mesh.', T_inv(T_sol_mesh, n_sol_mesh).')
% plot3(phi_sol_mesh(:,[1,end]).', n_sol_mesh(:,[1,end]).', ...
%     phi_sol_mesh(:,[1,end]).', 'r')
phi_approx_mesh = T_inv(T_mesh, n_mesh);
plot3(phi_mesh.', n_mesh.', phi_approx_mesh.')

% Either plot phi 0 to 1, or only valid range 0 to phi_max.
% plot3(phi_mesh(:,[1,end]).', n_mesh(:,[1,end]).', ...
%     phi_mesh(:,[1,end]).', 'r')

% Create linear indices to select maximum phi for each row in the mesh.
indices = sub2ind(size(phi_mesh), 1:size(phi_mesh, 1), ind.');
plot3([zeros(1, 100); phi_mesh(indices)], n_mesh(:,[1,end]).', ...
    [zeros(1, 100); phi_mesh(indices)], 'r')
plot3(phi_arr(ind), n_arr, phi_arr(ind), 'g-', 'LineWidth',2)


figure(4)
view(3)
hold on
xlabel("T")
ylabel("n")
zlabel("\phi_{pred} - \phi")

phi_ = T_inv(T_pred_mesh, n_pred_mesh);
res_phi = phi_ - phi_pred_mesh;
scatter3(T_pred_mesh, n_pred_mesh, ...
    reshape(output.residuals, size(T_pred_mesh)), "k.")
% residuals = output.residuals;
% plot( pressure,residuals,".")
% plot3(T_pred_mesh, n_pred_mesh, phi_sol)
% scatter3(phi_pred_mesh.', n_pred_mesh.', res_phi.', 'k.')
% plot3(phi_arr(ind), n_arr, zeros(size(n_arr)), 'g-', 'LineWidth',2)