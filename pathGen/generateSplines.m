clear; close all
load("referenceMCP.mat", "tMCP")
load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file

steps = 10;
tMCP_selection = tMCP(1:steps:end, :);
curve = cscvn(tMCP_selection.');


figure
hold on
fnplt(curve); 
plot(tMCP_selection(:,1), tMCP_selection(:,2), 'k.')
plot(t_r(:,1),t_r(:,2)) % plot right racetrack boundary
plot(t_l(:,1),t_l(:,2)) % plot left racetrack boundary

x_spline = ppmak(curve.breaks, curve.coefs(1:2:end, :), 1);
x_spline_ds = differentiate(x_spline);
x_spline_ds2 = differentiate(x_spline_ds);

y_spline = ppmak(curve.breaks, curve.coefs(2:2:end, :), 1);
y_spline_ds = differentiate(y_spline);
y_spline_ds2 = differentiate(y_spline_ds);

t_span = linspace(curve.breaks(1), curve.breaks(end), 100);
xi = ppval(x_spline,t_span);
yi = ppval(y_spline,t_span);
d1xi = ppval(x_spline_ds,t_span);
d1yi = ppval(y_spline_ds,t_span);
d2xi = ppval(x_spline_ds2,t_span);
d2yi = ppval(y_spline_ds2,t_span);
% K = abs(d1xi.*d2yi - d2xi.*d1yi) ./ sqrt(d1xi.^2+d1yi.^2).^3; 
K = (d1xi.*d2yi - d2xi.*d1yi) ./ sqrt(d1xi.^2+d1yi.^2).^3; % Signed curvature

tan_vec = [d1xi; d1yi] ./ vecnorm([d1xi; d1yi], 1);
ortho_vec = [cos(pi/2), -sin(pi/2); sin(pi/2), cos(pi/2)] * tan_vec;
ortho_vec_curve = K .* ortho_vec;

quiver(xi, yi, ortho_vec_curve(1,:), ortho_vec_curve(2,:), 1)
axis equal

function ppdf = differentiate(ppf)
% Spline Derivative
ppdf = ppf;
ppdf.order=ppdf.order-1;
ppdf.coefs=ppdf.coefs(:,1:end-1).*(ppdf.order:-1:1);
end