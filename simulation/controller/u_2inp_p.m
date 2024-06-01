function [U, debug] = u_2inp_p(X)
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
load("..\Die Controller ist sehr gut\referenceMCP.mat");
load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
t_r_x=t_r(:,1); % x coordinate of right racetrack boundary
t_r_y=t_r(:,2); % y coordinate of right racetrack boundary
t_l_x=t_l(:,1); % x coordinate of left racetrack boundary
t_l_y=t_l(:,2); % y coordinate of left racetrack boundary

% middle points of the line between coupled points
% t_m = 0.5 * (t_l + t_r);
t_m = tMCP;
dist = vecnorm([x, y] - t_m, 2, 2);
[~, pos_idx] = min(dist);
 % If there are multiple points, take the "next" one on the track
pos_idx = mod(pos_idx, 1950); % TODO
pos_next_idx = mod(pos_idx, 1950) + 1;
r = (t_l(pos_next_idx, :) - t_r(pos_next_idx, :)) ...
    / norm((t_l(pos_next_idx, :) - t_r(pos_next_idx, :)));
% tangent_vec = [0, 1; -1, 0] * r.';
% tanget_vec_forw = t_m(pos_next_idx, :) - t_m(pos_idx, :);

% psi_ref = atan2(tangent_vec(2), tangent_vec(1));
% 
% R_z_loc_glob = [cos(-pi/2+psi), -sin(-pi/2+psi); 
%                sin(-pi/2+psi), cos(-pi/2+psi)];
% p_loc_dot = R_z_loc_glob * [x_dot; y_dot];
% v_ref = 5;

x_loc = dot(r, [x, y] - t_m(pos_idx, :));
% x_loc_dot = p_loc_dot(1);
% 
% q_ref = [v_ref; beta_ref; ...
%     psi_ref; omega_ref; 0; 0];
% q = [v; beta; psi; omega; x_loc; x_loc_dot];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATE FEEDBACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("T_inv.mat", "T_inv");
persistent x_loc_err_int  x_loc_past
if isempty(x_loc_err_int)
    x_loc_err_int = 0;
    x_loc_past = 0;
end
x_loc_err = x_loc - 0;
x_loc_err_der = x_loc_err - x_loc_past;
x_loc_past = x_loc_err;

kp_v = 10000;
kp2_v = 10000;
kp_delta = 0.2;
ki_delta = -0.05;
kd_delta = 1;
u_delta = - kp_delta * x_loc_err - ki_delta * x_loc_err_int; %...
    % - kd_delta * x_loc_err_der;
delta = min(u_delta, 0.53);
delta = max(-0.53, delta);

x_loc_err_int = x_loc_err_int + 0.01 * x_loc_err;

if abs(delta) < 0.1
    v_ref = 10;
elseif abs(delta) < 0.15
    v_ref = 5;
elseif abs(delta) < 0.2
    v_ref = 3;
elseif abs(delta) < 0.4
    v_ref = 2;
else
    v_ref = 1.5;
end
v_err = v - v_ref;
u_v = - kp_v * v_err; % - kp2_v * x_loc_err;

if u_v >= 0
    R = 0.302;
    i0 = 3.91;
    i_g=[3.91 2.002 1.33 1 0.805];
    % T_want = (1/R * i0 * i_g(G))^-1 * u(1);
    [u1_max, G] = u1_max_f(v);
    u1_real = min(u_v, u1_max);
    T_real = (1/R * i0 * i_g(G))^-1 * u1_real;
    n = v * i_g(G) * i0 / R;
    phi = T_inv(T_real, n);
    Fb = 0;
else
    G = 1; % Does not matter
    phi = 0;
    Fb = -u_v;
end

zeta = 1;
% delta = 0;
% Fb = 0;
% phi = 1;

% delta=0; % steering angle
% G=1; %gear 1 ... 5
% Fb=0; % braking force
% zeta=0.5; %braking force distribution
% phi=0.2; % gas pedal position
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=[delta G Fb zeta phi]; % input vector
debug = [v_err x_loc_err x_loc_err_int x_loc_err_der];
end

