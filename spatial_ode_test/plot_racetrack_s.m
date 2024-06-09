%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot_racetrack
%
% plots the racetrack and your result
%
% files requested: racetrack.mat
%
% variables requested: Y (built by simulate_controlled_singletrack.m)
%
% The model is adopted from the scripts of M. Gerdts and D. Schramm,
% respectively.
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
f1 = figure('Name','racetrack','NumberTitle','off','Toolbar','figure','MenuBar','none');%,'OuterPosition',[0 0 460 1100]) % creates window for plot
hold on % allow for multiple plot commands within one figure
axis equal % eqal axis scaling
% axis([-50 70 -50 450]) % plot height and width
plot(t_r(:,1),t_r(:,2)) % plot right racetrack boundary
plot(t_l(:,1),t_l(:,2)) % plot left racetrack boundary
text(1,0,'\leftarrow finish/start','HorizontalAlignment','left') % finish/start annotation
plot(debug(:,1),debug(:,2),'r') % plot the x and y coordinates resulting fromy your controller
plot(debug(:,3),debug(:,4),'g')
xlabel('x') % label x axis
ylabel('y') % label y axies
% box % make a box around the plot

figure("Name", "States", "WindowState","maximized")
sgtitle("States")
subplot(6, 1, 1)
plot(sspan, Y(:,1))
ylabel("e_y")

subplot(6, 1, 2)
plot(sspan, Y(:,2))
ylabel("e_{\psi}")

subplot(6, 1, 3)
plot(sspan, Y(:,3))
ylabel("v")

subplot(6, 1, 4)
plot(sspan, Y(:,4))
ylabel("\beta")

subplot(6, 1, 5)
plot(sspan, Y(:,5))
ylabel("\omega")

subplot(6, 1, 6)
plot(sspan, Y(:,6))
ylabel("t")


% % U=[delta G Fb zeta phi];
figure("Name", "Controls", "WindowState","maximized")
tspan = sspan(1:end-1);
sgtitle("Control Inputs")
subplot(2, 1, 1)
plot(tspan, U(:, 1), "DisplayName", "u1")
ylabel("u1")

% subplot(4, 1, 2)
% plot(tspan, U(:, 2), "DisplayName", "F_b")
% ylabel("F_b")

% subplot(4, 1, 3)
% plot(tspan, U(:, 3), "DisplayName", "\zeta")
% ylabel("\zeta")

subplot(2, 1, 2)
plot(tspan, U(:, 2), "DisplayName", "\delta")
ylabel("\delta")

% 
% subplot(5, 1, 5)
% plot(tspan, U(:, 5), "DisplayName", "\phi")
% ylabel("\phi")
% xlabel("t")



% figure("Name", "Error", "WindowState","maximized")
% sgtitle("Errors")
% subplot(4,1,1)
% plot(tspan, debug(:, 1))
% ylabel("v_{err}")
% 
% subplot(4,1,2)
% plot(tspan, debug(:, 2))
% ylabel("x_{err}")
% 
% subplot(4,1,3)
% plot(tspan, debug(:, 3))
% ylabel("x_int_{err}")
% 
% subplot(4,1,4)
% plot(tspan, debug(:, 4))
% ylabel("x_der_{err}")

figure(f1)