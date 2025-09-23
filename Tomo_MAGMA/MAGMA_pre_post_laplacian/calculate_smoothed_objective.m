% calculate_smoothed_objective.m

function F_mu = calculate_smoothed_objective(x, L, lambda, params)
global N_vec
% Purpose:
%   Calculates the value of the smoothed objective function, F_mu(x).
%   F_mu(x) = f(x) + g_mu(x)

if isfield(params, 'mu')
    mu = params.mu;
else
    mu = 1e-6; % A reasonable default value
end

% Smooth part f(x)
f_val = (1/size(L,1))^2*0.5 * norm(L * x )^2;

% Smoothed L1-norm g_mu(x) from equation (2.3)
g_mu_val = lambda * sum(sqrt(mu^2 + x.^2));

F_mu = f_val + g_mu_val;

end