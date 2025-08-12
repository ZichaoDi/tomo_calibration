function [x, lambda, history] = solveBarrier(A, b, regWt, x0, mu, verbose)
%solver.solveBarrier() solve use Page 569, 11.3.1 Algorithm 11.1 the barrier method of the Byod@Stanford's bookd Convex Optimization.
%This function is only good for test for small cases where A is a matrix form and newton direction is solved eaxctly
% arg min_x = 0.5*||A*x - b||_2^2 + regWt * regFun(x),
%           = 0.5*||A*x - b||_2^2 + regWt * ||x||_1
%       s.t. x > 0


%
if nargin == 0
    [x, lambda, history] = testMe();
    return
end
if ~exist('mu', 'var') || isempty(mu)
    mu = 100;
end
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
%% note we combine [x; v] = x; and [rD; rP] = r. where x = x; rD = r(1:n);  v = x(n+1:end); rP = r(n+1:end)
[m, n]  = size(A);
tol = 1e-6;
t0 = 1;
% t0 = m/tol; % directly solve with one inner stop, will fail as (1) Matrix close to singular (2) Newton fails for non-smooth Hessian
objFun = @(x) 0.5*norm(A*x-b,2)^2 + regWt*sum(abs(x));  % 2018.02.05_added abs(x) for future to extend from non-negative x to real/complex x.;
history = [];
t = t0;
x = x0;
for kC = 1:100
    % centering step
    [x, retNt, k] = solver.solveNewtonExactly(A, b, regWt, x, t, verbose);

    history(1, kC) = k;
    history(2, kC) = m/t;
    history(3, kC) = objFun(x);
    fprintf('barrier: t=%.2e, dual gap = %.2e, Newton iteration=%d. objFunValwithT = %.2e, objFunVal = %.2e.\n', t, m/t, k, retNt.objFunVal(end), objFun(x));

    if m/t <= tol
        break
    end
    t = t*mu;
end

lambda = 1/t * 1./x;

if verbose
    [xx, yy] = stairs(cumsum(history(1,:)), history(2,:));
    [~, yy2] = stairs(cumsum(history(1,:)), history(3,:));
    figure;
    subplot(121); semilogy(xx, yy); title(sprintf('mu=%.2e, objFunVal = %.8e', mu, objFun(x)));xlabel('Newton iterations');  ylabel('duality gap');
    subplot(122); semilogy(xx, yy2-history(3,end)); xlabel('Newton iterations');  ylabel('objFunVal - objFunValOpt');
end
end

function [x, lambda, history] = testMe()
m = 100;
n = 500;
rng(0, 'twister');
A = rand(m, n);
rng(0, 'twister');
xfeasible = rand(n,1);
b = A*xfeasible;
rng(0, 'twister');
regWt = rand(1,1);

%x0 = [xfeasible; zeros(m, 1)];
%rng(0, 'twister');  x0 = [rand(n,1); ones(m, 1)];
%x0 = ones(n,1);
x0 = xfeasible;
%%
for mu = [2, 50, 150, 1000]
    [x, lambda, history] =  solveBarrier(A, b, regWt, x0, mu, true);
end
end
