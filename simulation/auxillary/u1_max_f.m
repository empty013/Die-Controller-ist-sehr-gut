function [u1_max, G] = u1_max_f(v_arr, G)
    fitting_coeffs = ...
    [2.80312112707768	-72.2104474231297	112.103949325811	-265.430812847204	40748.5346610902;
    4.54494664413907	-205.626417890814	3361.16520793424	-24485.4446199447	86634.5053493282;
    1.32459275849175	-94.5014221192759	2490.90200308857	-29236.3959382549	141776.206325086;
    0.814383300547809	-78.4260049045966	2814.15500610824	-44883.0445970889	278078.654116093;
    -0.00164595244592678	0.104047716958146	-6.70953768985322	108.649979111722	7599.44538850200];
    
    v_steps = [7.75757575757576	14.5454545454545	20.8484848484848	27.1515151515152	47.5151515151515];
    if exist('G','var')
     % third parameter does not exist, so default it to something
        u1_max = polyval(fitting_coeffs(G, :), v_arr);
        return
    end

    u1_max = zeros(size(v_arr));
    G = zeros(size(v_arr));
    for v_idx = 1:length(v_arr)
        v = v_arr(v_idx);
        if v <= v_steps(1)
            u1_max(v_idx) = polyval(fitting_coeffs(1, :), v);
            G(v_idx) = 1;
        elseif v <= v_steps(2)
            u1_max(v_idx) = polyval(fitting_coeffs(2, :), v);
            G(v_idx) = 2;
        elseif v <= v_steps(3)
            u1_max(v_idx) = polyval(fitting_coeffs(3, :), v);
            G(v_idx) = 3;
        elseif v <= v_steps(4)
            u1_max(v_idx) = polyval(fitting_coeffs(4, :), v);
            G(v_idx) = 4;
        else
            u1_max(v_idx) = polyval(fitting_coeffs(5, :), v);
            G(v_idx) = 5;
        end
    end
end