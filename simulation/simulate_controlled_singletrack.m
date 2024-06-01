function simulate_controlled_singletrack(t_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function simulate_controlled_singletrack(t_f)
%
% integrates the controlled single-track model until time t_f
%
% input: t_f (simulation time)
%
% files requested: racetrack.m ; singletrack.m ; ode1.m ; plot_racetrack.m
%
% plots built: racetrack
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
racetrack % builds the racetrack and saves it as racetrack.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_0=[-2.5;0;0;0;pi/2;0;0;0;0;0]; % initial value for integration
% X_0=[-2.5;240;0;0;pi/2;0;0;0;0;0]; % initial value for integration
% X_0=[22.5;340;0;0;3*pi/2+0.07;0;0;0;0;0]; % initial value for integration
tspan = 0:0.01:t_f;
[Y, U, debug] =ode1(@singletrack,tspan,X_0); % integrate with step zise 0.001

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_racetrack % plots the racetrack and your result
end