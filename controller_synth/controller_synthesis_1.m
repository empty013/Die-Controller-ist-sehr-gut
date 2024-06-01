clear 
syms x y v beta psi omega
syms x_dot y_dot
% syms m g l_r l_l r0 r1 r4
syms u1 u2 u3
syms delta zeta F_b
q = [x; y; v; beta; psi; x_dot; y_dot; omega];
u = [u1; F_b; zeta; delta];

m=1239; % vehicle mass
g=9.81; % gravitation
l_f=1.19016; % distance of the front wheel to the center of mass 
l_r=1.37484; % distance of the rear wheel to the center of mass
l = l_r + l_f;
R=0.302; % wheel radius
I_z=1752; % vehicle moment of inertia (yaw axis)
I_R=1.5; % wheel moment of inertia
r0=0.009; % coefficient (friction)
r1=0.002; % coefficient (friction)
r4=0.0003; % coefficient (friction)

B_f=10.96; % stiffnes factor (Pacejka) (front wheel)
C_f=1.3; % shape factor (Pacejka) (front wheel)
D_f=4560.4; % peak value (Pacejka) (front wheel)
E_f=-0.5; % curvature factor (Pacejka) (front wheel)
B_r=12.67; %stiffnes factor (Pacejka) (rear wheel)
C_r=1.3; %shape factor (Pacejka) (rear wheel)
D_r=3947.81; %peak value (Pacejka) (rear wheel)
E_r=-0.5; % curvature factor (Pacejka) (rear wheel)

dx_dt = v * cos(psi - beta);
dy_dt = v * sin(psi - beta);
psi_dot = omega;


mu = r0 + r1 * v + r4 * v^4;
a_f = delta - atan((l_f * omega - v * sin(beta)) / v * cos(beta));
a_r = atan((l_r * psi_dot + v * sin(beta)) / v * cos(beta));
F_x_f = - ((1 - zeta) * F_b + mu * m * g * l_r / l);
F_x_r = u1 * - (zeta * F_b + mu * m * g * l_r / l);
F_y_f = D_f * sin(C_f ...
    * atan(B_f *a_f - E_f * (B_f * a_f - atan(B_f * a_f))));
F_y_r = D_r * sin(C_r ...
    * atan(B_r *a_r - E_r * (B_r * a_r - atan(B_r * a_r))));


v_dot = 1/m * (F_x_r * cos(beta) + F_x_f * cos(delta + beta)...
            - F_y_r * sin(beta) - F_y_f * sin(delta + beta));
beta_dot = psi_dot - 1 / (m*v) ...
    * (F_x_r * sin(beta) + F_x_f * sin(delta + beta) ...
       + F_y_r * cos(beta) + F_y_f * cos(delta + beta));
psi_dot2 = 1/I_z * (F_y_f * l_f * cos(delta) - F_y_r * l_r ...
                    + F_x_f * l_f * sin(delta));
x_dot2 = 1/m * (...
    (F_x_r + F_x_f * cos(delta) - F_y_f * sin(delta)) * cos(psi) ...
    - (F_y_r + F_x_f * sin(delta) + F_y_f * cos(delta)) * sin(psi));
y_dot2 = 1/m * (...
    (F_x_r + F_x_f * cos(delta) - F_y_f * sin(delta)) * sin(psi) ...
    + (F_y_r + F_x_f * sin(delta) + F_y_f * cos(delta)) * cos(psi));

% q = [x; y; v; beta; psi; x_dot; y_dot; omega];
q_dot = [dx_dt; dy_dt; v_dot; beta_dot; psi_dot; x_dot2; y_dot2; psi_dot2];

syms x0 y0 psi0 x_dot0 y_dot0
q0 = [x0; y0; 5; 0; psi0; x_dot0; y_dot0; 0];
% u = [u1; F_b; zeta; delta];
u0 = [1; 0; 0.5; 0];

q_dot_jacobian_q = jacobian(q_dot, q);
A_sym = vpa(simplify(subs(q_dot_jacobian_q, [q; u], [q0; u0])));

q_dot_jacobian_u = jacobian(q_dot, u);
B_sym = vpa(simplify(subs(q_dot_jacobian_u, [q; u], [q0; u0])));

save("..\template_Matlab\system.mat", "A_sym", "B_sym")
% Kal = simplify([B, A*B, A^2*B, A^3*B, A^4*B, A^5*B, A^6*B, A^7*B]);
% rank(Kal);

n = size(A_sym, 1);
m = size(B_sym, 2);
Q = eye(n);
R = eye(m);
C = eye(n);
D = zeros(n, m);
psi0_num = pi/2;
A = double(subs(A_sym, psi0, psi0_num));
B = double(subs(B_sym, psi0, psi0_num));
sys = ss(A, B, C, D, "StateName", ...
    {"x", "y", "v", "beta", "psi", "x_dot", "y_dot", "omega"}, ...
    "InputName", {"u1", "Fbar", "zeta", "delta"});
[Abar,Bbar,Cbar,T,k] = ctrbf(A,B,C);
% [Abar,Bbar,Cbar,T,k] = ctrbf(sys);

% syms x_ref y_ref psi_ref
% v_ref = 5;
% x_dot_ref = v_ref * cos(psi_ref);
% y_dot_ref = v_ref * sin(psi_ref);
% beta_ref = 0;
% omega_ref = 0;
% q_ref = [x_ref; y_ref; v_ref; beta_ref; ...
%     psi_ref; omega_ref; x_dot_ref; y_dot_ref];
% [K, S, P] = lqr(A, B, Q, R);