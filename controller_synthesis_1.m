syms x y v beta psi omega
syms x_dot y_dot
% syms m g l_r l_l r0 r1 r4
syms u1 u2 u3
syms delta xi F_b
q = [x; y; v; beta; psi; omega; x_dot; y_dot];
u = [u1; F_b; xi; delta];

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

% x_dot = v * cos(psi - beta);
% y_dot = v * sin(psi - beta);
psi_dot = omega;


mu = r0 + r1 * v + r4 * v^4;
a_f = delta - atan((l_f * omega - v * sin(beta)) / v * cos(beta));
a_r = atan((l_r * psi_dot + v * sin(beta)) / v * cos(beta));
F_x_f = - ((1 - xi) * F_b + mu * m * g * l_r / l);
F_x_r = u1 * - (xi * F_b + mu * m * g * l_r / l);
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

q_dot = [x_dot; y_dot; v_dot; beta_dot; psi_dot; x_dot2; y_dot2; psi_dot2];

syms x0 y0 psi0 x_dot0 y_dot0
q0 = [x0; y0; 5; 0; psi0; 0; x_dot0; y_dot0];
u0 = [0; 0; 0.5; 0];

q_dot_jacobian_q = jacobian(q_dot, q);
A = simplify(subs(q_dot_jacobian_q, [q; u], [q0; u0]));

q_dot_jacobian_u = jacobian(q_dot, u);
B = simplify(subs(q_dot_jacobian_u, [q; u], [q0; u0]));
% K = simplify([B, A*B, A^2*B, A^3*B, A^4*B, A^5*B, A^6*B, A^7*B]);
% rank(K);

n = size(A, 1);
m = size(B, 2);
Q = eye(n);
R = eye(m);
[K, S, P] = lqr(A, B, Q, R);