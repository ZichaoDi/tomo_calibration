function F_mu = calculate_smoothed_objective(x, obj_f, lambda, mu)
% Calculates the value of the smoothed objective function
% F_mu(x) = f(x) + g_mu(x)
% where g_mu is the smoothed L1-norm
    
    f_val = obj_f(x);
    g_mu_val = lambda * sum(sqrt(mu^2 + x.^2));
    F_mu = f_val + g_mu_val;
end