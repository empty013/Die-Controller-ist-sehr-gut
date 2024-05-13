clear; close all
syms alpha_f
Bf =    10.96;              %stiffness factor
Cf =    1.3;                %shape factor
Df =    4560.4;             %peak value
Ef =    -0.5;               %curvature factor
figure(1)
hold on
alpha_f_arr = linspace(-pi/4, pi/4);

plot(alpha_f_arr, lateralFrontForce(alpha_f_arr))
plot(alpha_f_arr, lateralFrontForce_lin(alpha_f_arr))
grid on
xlabel("\alpha_f")
ylabel("F_{y,f}")
ylim([-5e3, 5e3])

F_yf = Df*sin(Cf*atan(Bf*alpha_f - Ef*(Bf*alpha_f - atan(Bf*alpha_f))));
F_yf_dot_sym = jacobian(F_yf);
F_yf_dot = matlabFunction(simplify(subs(F_yf_dot_sym, alpha_f, 0)) * alpha_f);
% plot(alpha_f_arr, F_yf_dot(alpha_f_arr))


function F_yf = lateralFrontForce(alpha_f)
    Bf =    10.96;              %stiffness factor
    Cf =    1.3;                %shape factor
    Df =    4560.4;             %peak value
    Ef =    -0.5;               %curvature factor
    % alpha_f = delta - arctan((lf*psi_dot - v*sin(beta))/(v*cos(beta)));
    F_yf = Df*sin(Cf*atan(Bf*alpha_f - Ef*(Bf*alpha_f - atan(Bf*alpha_f))));
end

function F_yf = lateralFrontForce_lin(alpha_f)
    Bf =    10.96;              %stiffness factor
    Cf =    1.3;                %shape factor
    Df =    4560.4;             %peak value
    Ef =    -0.5;               %curvature factor
    % alpha_f = delta - arctan((lf*psi_dot - v*sin(beta))/(v*cos(beta)));
    F_yf = Df*(Cf*(Bf*alpha_f - Ef*(Bf*alpha_f - (Bf*alpha_f))));
end