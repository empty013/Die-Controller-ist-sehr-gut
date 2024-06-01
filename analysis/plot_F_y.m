clear all
%%Parameters
m  =    1239;               %mass 
g  =    9.81;               %gravity constant
lf =    1.19016;            %front length
lr =    1.37484;            %rear length
R  =    0.302;              %wheel radius
Iz =    1752;               %vehicle moment of inertia
IR =    1.5;                %wheel moment of inertia 
i0 =    3.91;               %motor transmission
i1 =    3.91;               %first gear transmission
i2 =    2.002;              %second gear transmission
i3 =    1.33;               %third gear transmission
i4 =    1;                  %fourth gear transmission
i5 =    0.805;              %fifth gear transmission
i  =    [i0,i1,i2,i3,i4,i5];%gear transmissions
Bf =    10.96;              %stiffness factor
Cf =    1.3;                %shape factor
Df =    4560.4;             %peak value
Ef =    -0.5;               %curvature factor
Br =    12.67;              %stiffness factor
Cr =    1.3;                %shape factor
Dr =    3947.81;            %peak value
Er =    -0.5;               %curvature factor
p0 =    -37.8;              %coefficient
p1 =    1.54;               %coefficient
p2 =    -0.0019;            %coefficient
q0 =    -34.9;              %coefficient
q1 =    -0.04775;           %coefficient
r0 =    0.009;              %coefficient
r1 =    7.2*(10^-5);        %coefficient
r2 =    0;                  %coefficient
r3 =    0;                  %coefficient
r4 =    5.0388*(10^-10);    %coefficient


%%
% alpha_arr = linspace(-1,1)';
% plot(alpha_arr, lateralRearForce(alpha_arr))

v_arr = linspace(0, 20)';
plot(v_arr, m*g*lr/(lr+lf)*friction(v_arr))

function F_yr = lateralRearForce(alpha_r)
    Br =    12.67;              %stiffness factor
    Cr =    1.3;                %shape factor
    Dr =    3947.81;            %peak value
    Er =    -0.5;               %curvature factor
    F_yr = Dr*sin(Cr*atan(Br*alpha_r ...
        - Er*(Br*alpha_r - atan(Br*alpha_r))));
end

function mu = friction(v)
    r0 =    0.009;              %coefficient
    r1 =    7.2*(10^-5);        %coefficient
    r2 =    0;                  %coefficient
    r3 =    0;                  %coefficient
    r4 =    5.0388*(10^-10);    %coefficient
    mu = r0 + r1*abs(v) + r2*abs(v).^2 + r3*abs(v).^3 + r4*abs(v).^4;
end