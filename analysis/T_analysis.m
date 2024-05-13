clear; close all
v_min = 0;
v_max = 40;
v_res = 50;
phi_res = 50;
v_arr = linspace(v_min, v_max, v_res);
phi_arr = linspace(0, 1, phi_res);
[phi_mesh, v_mesh]  = meshgrid(phi_arr, v_arr);

syms v phi
R = 0.302;
i0 = 3.91;
i_g=[3.91 2.002 1.33 1 0.805];

figure(1)
view(3)
xlabel("\phi")
ylabel("v")
zlabel("F_x")
hold on

figure(2)
view(3)
xlabel("\phi")
ylabel("v")
zlabel("T_M")
hold on

G_max = 5;
colormap = hsv(G_max);
u_max = 0;
T_max = 0;
for G = 1:G_max
    n = v * i_g(G) * i0 / R;
    T_M=200*phi*(15-14*phi)- ...
        200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); 
    T_M_f = matlabFunction(T_M, "Vars", [phi, v]);
    T_mesh = T_M_f(phi_mesh, v_mesh);
    % u1 = 1/R * i_g(G) * i0 * T_M_f;
    % u1 = matlabFunction(u1_sym, "Vars", [phi, v]);
    u1_mesh = 1/R * i_g(G) * i0 * T_mesh;
    color = colormap(G,:);

    figure(1);
    h1 = surf(phi_mesh, v_mesh, u1_mesh, "FaceAlpha", "0.5");
    u_max = max(u_max, max(u1_mesh, [], "all"));
    set(h1, "FaceColor", color, "LineStyle","none");

    figure(2);
    h1 = surf(phi_mesh, v_mesh, T_mesh, "FaceAlpha", "0.5");
    T_max = max(T_max, max(T_mesh, [], "all"));
    set(h1, "FaceColor", color, "LineStyle","none");
end

figure(1)
zlim([0, u_max])
legend show

figure(2)
zlim([0, T_max])
legend show