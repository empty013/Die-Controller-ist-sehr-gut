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

%%Initial Values
x0 =        -2.5;   %longitudinal position
y0 =        0;      %lateral position
v0 =        0;      %velocity
beta0 =     0;      %side-slip angle
psi0 =      pi/2;   %yaw angle
omega0 =    0;      %yaw rate
x_dot0 =    0;      %longitudinal velocity
y_dot0 =    0;      %lateral velocity
psi_dot0 =  0;      %yaw rate (redundant)
phi_dot0 =  0;      %wheel rotary frequency

%%Input Constraints
delta_min = 0.53;   %min steering angle
delta_max = 0.53;   %max steering angle
Fb_min =    0;      %min breaking force
Fb_max =    15000;  %max breaking force
zeta_min =  0;      %min breaking distribution
zeta_max =  1;      %max breaking distribution
phi_min =   0;      %min gas pedal position
phi_max =   1;      %max gas pedal position
G = [1,2,3,4,5];    %gear positions

%%Functions
%inertial coordinates's rate of change
function [x_dot, y_dot] = xy_dot(v, psi, beta)
    x_dot = v*cos(psi-beta);
    y_dot = v*sin(psi-beta);
end

%velocity
function v = velocity(x_dot, y_dot)
    v = sqrt(pow(x_dot,2) + pow(y_dot,2));
end

%acceleration
function v_dot = acceleration(F_xf, F_xr, F_yf, F_yr, beta, delta)
    v_dot = 1/m*(F_xr*cos(beta) + F_xf*cos(beta+delta) - F_yr*sin(beta) - F_yf*sin(delta+beta));
end

%friction
function mu = friction(v)
    mu = r0 + r1*abs(v) + r2*pow(abs(v),2) + r3*pow(abs(v),3) + r4*pow(abs(v),4);
end

%longitudinal rear force
function F_xr = longitudinalRearForce(Fb, zeta, beta, v)
    signum = sign(v*cos(beta));
    F_br = -Fb*(1-zeta);                %breaking force rear
    F_mur = -friction(v)*(m*g*lr/l);    %friction force rear
    F_xr = signum*(F_br + F_mur);
end

%longitudinal front force
function F_xf = longitudinalFrontForce(phi, Fb, zeta, k, beta, v)
    signum = sign(v*cos(beta));
    F_bf = -Fb*zeta;                        %breaking force front
    F_muf = -friction(v)*(m*g*lf/l);        %friction force front
    F_M  = 1/R*i(G(k))*i(0)*T_M(phi, G(k)); %motor force
    F_xf = signum*(F_M + F_bf + F_muf);
end

%lateral rear force
function F_yr = lateralRearForce()
    alpha_r = atan((lr*psi_dot - v*sin(beta))/(v*cos(beta)));
    F_yr = Dr*sin(Cr*atan(Br*alpha_r - Er(Br*alpha_r - atan(Br*alpha_r))))
end

%lateral front force
function F_yf = lateralFrontForce(psi_dot, v, beta)
    alpha_f = delta - atan((lf*psi_dot - v*sin(beta))/(v*cos(beta)));
    F_yf = Df*sin(Cf*atan(Bf*alpha_f - Ef*(Bf*alpha_f - atan(Bf*alpha_f))));
end

%motor torque
function T_M = motorTorque(phi, k, v)
    n = motorRotaryFrequency(k, v);
    T_M = 200*phi*(15-14*phi)*(1 - pow((30/pi*n/4800),5*phi))
end

%motor rotary frequency
function n = motorRotaryFrequency(k, v)
    n = v*i(G(k))*i0/R;
end
