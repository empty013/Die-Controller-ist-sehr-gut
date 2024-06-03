function [k_ref, psi_ref, x_ref, y_ref] = referencePath(s)
load("referenceSplines.mat", "x_spline", "x_spline_ds", "x_spline_ds2", ...
    "y_spline", "y_spline_ds", "y_spline_ds2");
s = s + 100; % <-------------------------------------------------REMOVE-------------<
x_ref = ppval(x_spline,s);
y_ref = ppval(y_spline,s);
x_dot_ref = ppval(x_spline_ds,s);
y_dot_ref = ppval(y_spline_ds,s);
x_dot2_ref = ppval(x_spline_ds2,s);
y_dot2_ref = ppval(y_spline_ds2,s);
k_ref = (x_dot_ref.*y_dot2_ref - x_dot2_ref.*y_dot_ref) ...
    ./ sqrt(x_dot_ref.^2+y_dot_ref.^2).^3; % Signed curvature

% n = [x_dot_ref; y_dot_ref] / norm([x_dot_ref; y_dot_ref]);
psi_ref = atan2(y_dot_ref, x_dot_ref);
end

