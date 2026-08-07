function [f_val, grad_val] = gradient_checker_wrapper(x, obj_func, grad_func)
%GRADIENT_CHECKER_WRAPPER Helper function to format objective and gradient for checkGradients.
    f_val = obj_func(x);
    grad_val = grad_func(x);
end