function s = find_s(x0)
load("referenceSplines.mat", "L")
s_grid = linspace(0, L, 100000);
[~, ~, x_arr, y_arr] = referencePath(s_grid);
p_arr = [x_arr; y_arr];
diff_arr = vecnorm(p_arr - x0, 2, 1);
[~, I] = min(diff_arr);
s = s_grid(I);
end