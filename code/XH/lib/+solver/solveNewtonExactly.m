function [x, ret, k] = solveNewtonExactly(AMatrix, b, regWt, x0, t, verbose)
%solver.solveNewtonExactly(): solve using page 487, 9.5.2 Newton's method of the Byod@Stanford's bookd Convex Optimization.
%This solve Newton direction exactly, which is only good for small test cases, where AMatrix is a full matrix
%  arg min_x = 0.5*||A*x - b||_2^2 + regWt * regFun(x) + 1/t * Phi(x), 
%           = 0.5*||A*x - b||_2^2 + regWt * ||x||_1 + 1/t * (-sum(log x)), 
% 
if nargin == 0
    [x, ret, k] = testMe();
    return;
end

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = true;
end

if ~exist('t', 'var') || isempty(t)
    t = 1;
end

%% 
%tol = max(1e-14*t, 1e-8);
tol = 1e-8;
[m, n] = size(AMatrix);
ATA = AMatrix.'*AMatrix;


gradFun = @(x) AMatrix.'*(AMatrix*x-b) + regWt*ones(n,1) -  1/t * 1./x;
hessianFun = @(x) ATA + 1/t * diag(1./(x.^2));
objFunT = @(x) 0.5*norm(AMatrix*x-b, 2)^2 + regWt*sum(x) - 1/t * sum(log(x));
objFun = @(x) 0.5*norm(AMatrix*x-b,2)^2 + regWt*sum(abs(x));  % 2018.02.05_added abs(x) for future to extend from non-negative x to real/complex x.;

theta = 0.01;
zeta = 0.5;

k=0; 
x = x0;
ret = [];
ret.objFunValT(k+1) = objFunT(x);
ret.objFunVal(k+1) = objFun(x);
for k = 1:1000    
    % 1. compute the primal and dual Newton steps dx, dv   
    % Method I: directly solve the inverse     
    grad = gradFun(x);
    hessian = hessianFun(x);    
    dx = -hessian\grad;
    
	% Method II: solve by block elimination, note the Hessian is a diagonal matrix, so easily invertable
%     HinvDiag = x.^2;
%     dx = zeros(n+m, 1, 'like', x);
%     dx(n+1:end) = (AMatrix*bsxfun(@times, HinvDiag, AMatrix.')) \ (r(n+1:end) - AMatrix*(HinvDiag.*r(1:n)));
%     dx = - HinvDiag .*  (AMatrix.'*dx(n+1:end) + r(1:n));
    
    ret.ntDecrementSquare(k) = - grad'*dx;
    if verbose
        fprintf('Iteration %d, objFunVal = %15.8e, objFunValT = %.4e, Newton decrement square = %.4e.\n', k-1, objFun(x), ret.objFunVal(k), ret.ntDecrementSquare(k));
    end
    
    % Converged to feasible solution && small residual
    if ret.ntDecrementSquare(k)/2 <= tol
        break
    end

    % 2. Backtracking line search    
    Adx = AMatrix*dx;
    r = b - AMatrix*x;    
    objFunT_x2s = @(x2, alpha) (0.5*norm(r(:)-alpha*Adx(:),2)^2 + regWt*sum(x2(:))) - 1/t*sum(log(x2(:)));    % x2 = x+alpha*dx, r2 = b-AMatrix*x2 = r - alpha*Adx
    [~, x2, objFunValTx2] = backTrackLineSearch(x, dx, -ret.ntDecrementSquare(k), objFunT_x2s,  ret.objFunValT(k), theta, zeta, verbose); 
    
    % 3. Update x and corresponding function values
    x = x2;    
    ret.objFunValT(k+1) = objFunValTx2; % = objFunT(x);
    ret.objFunVal(k+1) = objFun(x);
end


end



function [x, ret, k] = testMe()
%%
m = 100;
n = 500;
rng(0, 'twister');
AMatrix = rand(m, n);
rng(0, 'twister');
xfeasible = rand(n,1);
b = AMatrix*xfeasible;
rng(0, 'twister');
regWt = rand(1,1);

%x0 = [xfeasible; zeros(m, 1)];
%rng(0, 'twister');  x0 = [rand(n,1); ones(m, 1)];
%x0 = ones(n,1); % ones(n,1) >  0.001*ones(n,1) >  xfeasible >> eps*ones(n,1), where ">" means better, ">>" means far better
x0 = xfeasible;
t = 1e-1; 
verbose = true;
[x, ret, k] = solver.solveNewtonExactly(AMatrix, b, regWt, x0, t, verbose);
figure; %multAxes(@semilogy, {ret.ntDecrementSquare, ret.objFunValT-ret.objFunValT(end), ret.objFunVal-ret.objFunVal(end)})
subplot(131); semilogy(ret.ntDecrementSquare); title('Newton Decreament Square');
subplot(132); semilogy(ret.objFunValT-ret.objFunValT(end)); title('objeFunValT-objFunValTBest');
subplot(133); semilogy(ret.objFunVal-ret.objFunVal(end)); title('objeFunVal-objFunValBest');
end
