function [U] = controller(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function [U] = controller(X)
%
% controller for the single-track model
%
% inputs: x (x position), y (y position), v (velocity), beta
% (side slip angle), psi (yaw angle), omega (yaw rate), x_dot (longitudinal
% velocity), y_dot (lateral velocity), psi_dot (yaw rate (redundant)), 
% varphi_dot (wheel rotary frequency)
%
% external inputs (from 'racetrack.mat'): t_r_x (x coordinate of right 
% racetrack boundary), t_r_y (y coordinate of right racetrack boundary),
% t_l_x (x coordinate of left racetrack boundary), t_l_y (y coordinate of
% left racetrack boundary)
%
% outputs: delta (steering angle ), G (gear 1 ... 5), F_b (braking
% force), zeta (braking force distribution), phi (gas pedal position)
%
% files requested: racetrack.mat
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% prepared by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de
%
% written by *STUDENT*, *DATE*
% mailto:*MAILADDRESS*


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% state vector
x=X(1); % x position
y=X(2); % y position
v=X(3); % velocity (strictly positive)
beta=X(4); % side slip angle
psi=X(5); % yaw angle
omega=X(6); % yaw rate
x_dot=X(7); % longitudinal velocity
y_dot=X(8); % lateral velocity
psi_dot=X(9); % yaw rate (redundant)
varphi_dot=X(10); % wheel rotary frequency (strictly positive)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% racetrack
load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
t_r_x=t_r(:,1); % x coordinate of right racetrack boundary
t_r_y=t_r(:,2); % y coordinate of right racetrack boundary
t_l_x=t_l(:,1); % x coordinate of left racetrack boundary
t_l_y=t_l(:,2); % y coordinate of left racetrack boundary

% middle points of the line between coupled points
t_m = 0.5 * (t_l + t_r);
dist = vecnorm([x, y] - t_m, 2, 2);
[~, pos_idx] = min(dist);
 % If there are multiple points, take the "next" one on the track
pos_idx = pos_idx(end);
pos_next_idx = mod(pos_idx, 1950) + 1;
r = 0.5 * (t_l(pos_next_idx, :) - t_r(pos_next_idx, :));
tangent_vec = [0, 1; -1, 0] * r.';

psi_ref = atan2(tangent_vec(2), tangent_vec(1));

% x_ref = x;
% y_ref = y;
v_ref = 5;
beta_ref = 0;
omega_ref = 0;
x_dot_ref = v_ref * cos(psi_ref);
y_dot_ref = v_ref * sin(psi_ref);
q_ref = [v_ref; beta_ref; ...
    psi_ref; x_dot_ref; y_dot_ref; omega_ref];
q = [v; beta; psi; x_dot; y_dot; omega];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATE FEEDBACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms psi0
load("system.mat", "A_sym", "B_sym");
load("T_inv.mat", "T_inv");
n = size(A_sym, 1);
m = size(B_sym, 2);
Q = eye(n);
R = eye(m);
A = double(subs(A_sym, psi0, psi));
B = double(subs(B_sym, psi0, psi));
[K, ~, ~] = lqr(A, B, Q, R);

u = - K * (q - q_ref);
% u = [u1; F_b; zeta; delta];

R = 0.302;
i0 = 3.91;
i_g=[3.91 2.002 1.33 1 0.805];
G = 2;
n = v * i_g(G) * i0 / R;
phi = T_inv((1/R * i0 * i_g(G))^-1 * u(1), n);

Fb = u(2);
zeta = u(3);
delta = u(4);

% delta=0; % steering angle
% G=1; %gear 1 ... 5
% Fb=0; % braking force
% zeta=0.5; %braking force distribution
% phi=0.2; % gas pedal position
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=[delta G Fb zeta phi]; % input vector
end

