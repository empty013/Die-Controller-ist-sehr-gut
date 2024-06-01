clear; %close all
figure(1)
hold on

syms phi v
R = 0.302;
i0 = 3.91;
i_g=[3.91 2.002 1.33 1 0.805];

v_max = 48;
v_res = 100;
phi_res = 100;
v_arr = linspace(0, v_max, v_res);
phi_arr = linspace(0, 1, phi_res);
[phi_mesh, v_mesh]  = meshgrid(phi_arr, v_arr);

u1_mesh_arr = zeros([size(phi_mesh), 5]);

for G = 1:5
    n = v * i_g(G) * i0 / R;
    T_M=200*phi*(15-14*phi)- ...
        200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); 
    T_M_f = matlabFunction(T_M, "Vars", [phi, v]);
    u1_mesh = 1/R * i0 * i_g(G) * T_M_f(phi_mesh, v_mesh);
    u1_mesh_arr(:,:,G) = u1_mesh;
end
[u1_max_mesh, G_ind] = max(u1_mesh_arr, [], 3);
[u1_max_arr, phi_max_ind] = max(u1_max_mesh, [], 2);

% Create linear indices to select maximum G for each row in the mesh.
G_indices = sub2ind(size(G_ind), 1:size(G_ind, 1), phi_max_ind.');
G_arr = G_ind(G_indices);
G_step_idx = find(diff(G_arr));
G_step_idx = [G_step_idx, v_res-1];
G_start_idx = 1;
fitting_coeffs = zeros(length(G_step_idx), 5);
for G_idx = 1:length(G_step_idx)
    G_end_idx = G_step_idx(G_idx) + 1;
    v_data = v_arr(G_start_idx:G_end_idx);
    u1_data = u1_max_arr(G_start_idx:G_end_idx);
    phi_data = phi_arr(phi_max_ind(G_start_idx:G_end_idx));
    % plot3(phi_data, v_data, u1_data);
    G_start_idx = G_end_idx;

    [xData, yData] = prepareCurveData( v_data, u1_data );

    % Set up fittype and options.
    ft = fittype( 'poly4' );
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    fitting_coeffs(G_idx, :) = coeffvalues(fitresult);
    % plot3(phi_data, v_data, fitresult(v_data));
end
v_step = v_arr(G_step_idx);
save("u1_max_data", "v_step", "fitting_coeffs");

surf(phi_mesh, v_mesh, u1_max_mesh, "LineStyle","none");
% plot3(phi_arr(phi_max_ind), v_arr, u1_max_arr, 'r');
tmp = u1_max_f(v_arr);
plot3(phi_arr(phi_max_ind), v_arr, u1_max_f(v_arr));

plot3(zeros())

