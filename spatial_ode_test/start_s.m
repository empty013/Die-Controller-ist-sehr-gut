clear all; close all
addpath("..\pathGen\")
global X0 s_delta
% X = [e_y; e_psi; v; beta; omega; t];
[~, psi0_ref, x0_ref, y0_ref] = referencePath(0);
x0 = x0_ref;
y0 = y0_ref;
psi0 = psi0_ref;
% x0 = -2.5;
% y0 = 0;
% psi0 = pi/2;
% psi_ref = pi/2;
% x_ref = -2.5;
% y_ref = 0;
e_y_0 = cos(psi0_ref) * (y0 - y0_ref) - sin(psi0_ref) * (x0 - x0_ref);
e_psi_0 = psi0 - psi0_ref;
X0 = [e_y_0; e_psi_0; 10; 0; 0; 0];
s_delta = 0.01;
simulate_controlled_singletrack_s(10)