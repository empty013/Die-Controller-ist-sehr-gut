clear all; close all
addpath("controller\")
addpath("auxillary\")
addpath("..\pathGen\")
addpath("..\controller_synth\")
global controller_func X_0 t_delta
s0 = 0;
[~, psi0_ref, x0_ref, y0_ref] = referencePath(s0);
x0 = x0_ref-0.1;
y0 = y0_ref;
psi0 = psi0_ref;
v0 = 1;
X_0=[x0;y0;v0;0;psi0;0;0;0;0;0];
% X_0=[-2.5;0;1;0;pi/2;0;0;0;0;0];
% X_0=[-2.5;240;0;0;pi/2;0;0;0;0;0];
% X_0=[22.5;340;0;0;3*pi/2+0.07;0;0;0;0;0];
% % X_0=[-2.5;0;0;0;pi/2;0;0;0;0;0]; % initial value for integration

% controller_func = @u_2inp_p;
controller_func = @u_NMPC;
t_delta = 0.01;
tic
% simulate_controlled_singletrack(97.9)
simulate_controlled_singletrack(0.5)
toc