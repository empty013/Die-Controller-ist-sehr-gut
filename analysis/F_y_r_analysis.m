clear; close all
alpha_r_arr = linspace(-pi/4, pi/4);

plot(alpha_r_arr, lateralRearForce(alpha_r_arr))
grid on
xlabel("\alpha_r")
ylabel("F_{y,r}")


function F_yr = lateralRearForce(alpha_r)
    Br =    12.67;              %stiffness factor
    Cr =    1.3;                %shape factor
    Dr =    3947.81;            %peak value
    Er =    -0.5;               %curvature factor
    % alpha_r = arctan((lr*psi_dot - v*sin(beta))/(v*cos(beta)));
    F_yr = Dr*sin(Cr*atan(Br*alpha_r - Er*(Br*alpha_r - atan(Br*alpha_r))));
end