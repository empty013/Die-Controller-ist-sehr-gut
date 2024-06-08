% function u = synth_NMPC(s, q0)
clear; 
addpath("..\pathGen\")
% load("constraints", "ref_func", "x0_f", "lb_f", "ub_f");
load("constraints", "k_ref_func", "x0_f", "lb", "ub", "s_end", "s_delta");
opts = optimoptions('fmincon','Algorithm','interior-point', 'Display','none', ...
    "MaxFunctionEvaluations",3000, "UseParallel",true);%, 'SpecifyConstraintGradient',true);
% s_end = 0.5;
% s_delta = 0.01;
N = s_end / s_delta;

% s = 113;
s = 0;
% q0 = [0.105; 0.437; 3; 0; 0; 0];
q0 = [0.1; 0; 3; 0; 0; 0];
e_y = [q0(1)];
e_psi = [q0(2)];
v = [q0(3)];
beta = [q0(4)];
omega = [q0(5)];
t = [q0(6)];

u1 = [];
F_b = [];
zeta = [];
delta = [];
x_ref_arr = [];
y_ref_arr = [];
psi_ref_arr = [];
k_ref_arr = [];
sspan_arr = [];

iters = 5;
for ii = 1:iters
    ii
    k_ref = k_ref_func(s);
    x0 = x0_f(q0);
    % lb = lb_f(q0);
    % ub = ub_f(q0);
    f = @costFunction;
    nonlin = @(x) c_nonlin(x, q0, k_ref);
    % f = @(x) costFunction(x, ref);
    
    tic
    [sol, fval, ~, output] = fmincon(f, x0, [], [], [], [], lb, ub, nonlin, opts);
    toc
    
    start_state = 3;
    start_input = 1;
    stride = 8;
    e_y = [e_y; sol(start_state:stride:end)];
    e_psi = [e_psi; sol(start_state+1:stride:end)];
    v = [v; sol(start_state+2:stride:end)];
    beta = [beta; sol(start_state+3:stride:end)];
    omega = [omega; sol(start_state+4:stride:end)];
    t = [t; sol(start_state+5:stride:end)];
    
    u1 =  [u1; sol(start_input:stride:end)];
    % F_b = [F_b; sol(start_input+1:stride:end)];
    % zeta = [zeta; sol(start_input+2:stride:end)];
    delta = [delta; sol(start_input+1:stride:end)];
    
    sspan = linspace(s, s + s_end, N);
    sspan_arr = [sspan_arr(1:end-1), sspan];
    [k_ref, psi_ref, x_ref, y_ref] = referencePath(sspan);
    x_ref_arr = [x_ref_arr(1:end-1), x_ref];
    y_ref_arr = [y_ref_arr(1:end-1), y_ref];
    psi_ref_arr = [psi_ref_arr(1:end-1); psi_ref.'];
    k_ref_arr = [k_ref_arr(1:end-1); k_ref.'];
    
    % s = s + s_end;
    x_glob = x_ref_arr(end) - e_y(end) .* sin(psi_ref(end));
    y_glob = y_ref_arr(end) + e_y(end) .* cos(psi_ref(end));
    s = find_s([x_glob; y_glob]);
    
    q0 = [e_y(end); e_psi(end); v(end); beta(end); omega(end); t(end)];
    q0 = [e_y(end); e_psi(end); v(end); beta(end); omega(end); 0];
end
%%
close all
sspan = sspan_arr;
U_arr = [u1, F_b, zeta, delta];
save("input_recording", "U_arr")

x = x_ref_arr.' - e_y .* sin(psi_ref_arr);
y = y_ref_arr.' + e_y .* cos(psi_ref_arr);
% TODO: CHECK e_psi!!!
% Not matching with graph

load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
f1 = figure('Name','racetrack','NumberTitle','off','Toolbar','figure','MenuBar','none');%,'OuterPosition',[0 0 460 1100]) % creates window for plot
hold on
axis equal
plot(t_r(:,1),t_r(:,2)) % plot right racetrack boundary
plot(t_l(:,1),t_l(:,2)) % plot left racetrack boundary
plot(x,y,'r') % plot the x and y coordinates resulting fromy your controller
xlabel('x') % label x axis
ylabel('y') % label y axies
plot(x_ref_arr, y_ref_arr, 'g')

figure("Name", "States", "WindowState","maximized")
sgtitle("States")
subplot(6, 1, 1)
plot(sspan, e_y)
ylabel("e_y")

subplot(6, 1, 2)
plot(sspan, e_psi)
ylabel("e_{\psi}")

subplot(6, 1, 3)
plot(sspan, v)
ylabel("v")

subplot(6, 1, 4)
plot(sspan, beta)
ylabel("\beta")

subplot(6, 1, 5)
plot(sspan, omega)
ylabel("\omega")

subplot(6, 1, 6)
plot(sspan, t)
ylabel("t")



figure("Name", "Controls", "WindowState","maximized")

tspan = sspan(1:end-1);
sgtitle("Control Inputs")
subplot(2, 1, 1)
plot(tspan, u1, "DisplayName", "u1")
ylabel("u1")

% subplot(4, 1, 2)
% plot(tspan, F_b, "DisplayName", "F_b")
% ylabel("F_b")

% subplot(4, 1, 3)
% plot(tspan, zeta, "DisplayName", "\zeta")
% ylabel("\zeta")

subplot(2, 1, 2)
plot(tspan, delta, "DisplayName", "\delta")
ylabel("\delta")

figure(f1)