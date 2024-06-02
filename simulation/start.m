clear all; close all
addpath("controller\")
global controller_func X_0 t_delta
X_0=[-2.5;0;0;0;pi/2;0;0;0;0;0]; % initial value for integration
% X_0=[-2.5;240;0;0;pi/2;0;0;0;0;0];
% X_0=[22.5;340;0;0;3*pi/2+0.07;0;0;0;0;0];
controller_func = @u_2inp_p;
t_delta = 0.01;
simulate_controlled_singletrack(10)