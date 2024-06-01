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


% couple left and right track points
sizeL = size(t_l);
sizeR = size(t_r);
nL = sizeL(1,1);
nR = sizeR(1,1);
coupling = zeros(1,3);      % [indexLeftPoint, indexRightPoint, distance]
couplings = zeros(1950,3);  % [indexLeftPoint, indexRightPoint, distance]

for iL = 1:nL
    coupling(1,1) = iL;
    coupling(1,3) = power(10,6);
    pointL = t_r(iL,:);
    xL = pointL(1,1);
    yL = pointL(1,2);
    for iR = 1:nR
        pointR = t_l(iR,:);
        xR = pointR(1,1);
        yR = pointR(1,2);
        distance = power(xL-xR,2) + power(yL-yR,2);
        if (coupling(1,3) > distance && distance ~= 0)
            coupling(1,2) = iR;
            coupling(1,3) = distance; 
        end
    couplings(iL,:) = coupling;
    end
end


% middle points of the line between coupled points
mPoints = zeros(1950,4);
for iL = 1:nL
    iR = couplings(iL,2);
    mPoints(iL,1:2) = t_l(iL,:) - 0.5.*(t_r(iR,:) - t_l(iL,:));
    RotMatrix = [0 -1; 1 0];
    mPoints(iL,3:4) = 1/(couplings(iL,3)).*(t_r(iR,:) - t_l(iL,:))*RotMatrix;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATE FEEDBACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=0; % steering angle
G=1; %gear 1 ... 5
Fb=0; % braking force
zeta=0.5; %braking force distribution
phi=0.2; % gas pedal position
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=[delta G Fb zeta phi]; % input vector
end

