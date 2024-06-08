% function [c_ineq, c_eq] = c_nonlin(in1,in2, in3)
function [c_ineq, c_eq, c_ineq_grad, c_eq_grad] = c_nonlin(in1,in2,in3)
    c_ineq = [];
    c_eq = c_eq_f(in1, in2, in3);

    if nargout > 2
        c_ineq_grad = [];
        c_eq_grad = c_eq_dx_f(in1, in2, in3).';
    end
end
