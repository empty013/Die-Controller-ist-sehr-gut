clear;
syms beta v omega
syms u1 F_b zeta delta
syms x_ref y_ref k_ref psi_ref
syms e_y e_psi t
q = [e_y; e_psi; v; beta; omega; t];
% u = [u1; F_b; zeta; delta];
u = [u1; delta];
zeta = 1;
ref = [k_ref; psi_ref];
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


% psi = psi_ref + e_psi;
% psi_dot = omega;
% x_dot = v * cos(psi - beta);
% y_dot = v * sin(psi - beta);
x_dot = v * cos(beta);
y_dot = -v * sin(beta);

mu = r0 + r1 * v + r4 * v^4;
a_f = delta - atan((l_f * omega - v * sin(beta)) / (v * cos(beta)));
a_r = atan((l_r * omega + v * sin(beta)) / (v * cos(beta)));
F_x_f = -(1 - zeta) * 15000*F_b - mu * m * g * l_r / l;
% F_x_r = u1*10000 - zeta * 15000*F_b - mu * m * g * l_f / l;
F_x_r = 10000*u1 - mu * m * g * l_f / l;
F_y_f = D_f * sin(C_f ...
    * atan(B_f *a_f - E_f * (B_f * a_f - atan(B_f * a_f))));
F_y_r = D_r * sin(C_r ...
    * atan(B_r *a_r - E_r * (B_r * a_r - atan(B_r * a_r))));

% Linearizing
% q0 = [0; 0; 3; 0; 0; 0];
% u0 = zeros(size(u));
% F_y_f = subs(F_y_f, [q; u], [q0; u0]) ...
%     + subs(jacobian(F_y_f, [q; u]), [q; u], [q0; u0]) * [q; u];
% F_y_r = subs(F_y_r, [q; u], [q0; u0])...
%     + subs(jacobian(F_y_r, [q; u]), [q; u], [q0; u0]) * [q; u];


v_dot = 1/m * (F_x_r * cos(beta) + F_x_f * cos(delta + beta)...
            - F_y_r * sin(beta) - F_y_f * sin(delta + beta));
beta_dot = omega - 1 / (m*v) ...
    * (F_x_r * sin(beta) + F_x_f * sin(delta + beta) ...
       + F_y_r * cos(beta) + F_y_f * cos(delta + beta));
psi_dot2 = 1/I_z * (F_y_f * l_f * cos(delta) - F_y_r * l_r ...
                    + F_x_f * l_f * sin(delta));


% e_y = cos(psi_ref) * (y - y_ref) - sin(psi_ref) * (x - x_ref);
% phi_ref = atan2(y_ref, x_ref);
% e_phi = phi - phi_ref;
s_dot = 1 / (1 - e_y * k_ref) * (x_dot * cos(e_psi) - y_dot * sin(e_psi));
beta_dot = simplify(beta_dot / s_dot);
v_dot = simplify(v_dot / s_dot);
psi_dot2 = simplify(psi_dot2 / s_dot);
% psi_dot = psi_dot / s_dot;


e_psi_dot = omega / s_dot - k_ref;
e_y_dot = (x_dot * sin(e_psi) + y_dot * cos(e_psi)) / s_dot;
t_dot = 1 / s_dot;

f_sym = simplify([e_y_dot; e_psi_dot; v_dot; beta_dot; psi_dot2; t_dot]);

% Linearizing
% q0 = [0; 0; 3; 0; 0; 0];
% u0 = [0; 0; 0];
% f0 =  subs(f_sym, [q; u],  [q0; u0]);
% f_sym_lin = subs(jacobian(f_sym, [q; u]), [q; u], [q0; u0]);
% f_sym = f0 + f_sym_lin * [q; u];

save("f_s_sym", "f_sym")
f = matlabFunction(f_sym, "Vars", {q; u; ref}, "File","f_s");
% f = matlabFunction(f_sym, "Vars", {q; u; ref});
s_dot_f = matlabFunction(s_dot, "Vars", {q; u; ref});

% q0 = [0; 0; 1; 0; 0; 0];
% u0 = zeros(size(u));
% ref0 = [0; pi/2];

% f(q0, u0, ref0)
% s_dot_f(q0, u0, ref0)