function [x, objFunVals, cpuTimes, debias_start, SNRs] = solveDebias(b, AFun, AtFun, x, xGT, tolD, maxIterD, minIterD, isNonNegative, verbose, k, objFunVals, cpuTimes, SNRs, regWt, regFun, t0)
%--------------------------------------------------------------
% remove the bias from the l1 penalty, by applying CG to the least-squares problem 
% obtained by omitting the l1 term and fixing the zero coefficients at zero.
% Xiang changed Debias step as a seperate function, so that it could be used by TwIST or FISTA, as a second step of main
% algorithm
% TODO: Add non negative constraints here
% ============ Extra Step of Debias ==============================
% Note: Debiasing is an operation aimed at the computing the solution of the LS problem 
%                         arg min_x = 0.5*|| b - As * xs ||_2^2 
%                 where As is the  submatrix of A obatained by  deleting the columns of A corresponding of components
%                 of x set to zero by the TwIST algorithm                 
%  Compare our
%                         arg min_x = 0.5 xs'*As'*As*xs - xs'*As'*b
%  with that in https://en.wikipedia.org/wiki/Conjugate_gradient_method
%                         arg min_x = 0.5 x'*A*x - x'*b
% r = - grad(f(x)) = As'(b - As*xs) = As'*resid, where resid = b - As*xs
%  
%  'tolD' = stopping threshold/tolerance for the debiasing phase:
%
%  'maxIterD' = maximum number of iterations allowed in the debising phase of the algorithm.
%
%  'minIterD' = minimum number of iterations to perform in the debiasing phase of the algorithm.
%   
%   'x' = intial (sparse) estimate from Twist/Fista etc.
%   
%   'AFun':   b(:) = A*x(:)
%   'AtFun':  x(:) = A'*x(:)
%


solver.printAlgorithmBegin('algorithmName', mfilename);

idxZero =  abs(x) < 10*eps(class(x));
cont_debias_cg = 1;
debias_start = k;

% calculate initial residual
resid = b - AFun(x); %resid = AFun(x) - b;
r = AtFun(resid);

% mask out the zeros
r(idxZero) = 0;
rTr = r(:)'*r(:);

% set convergence threshold for the residual || RW x - b ||_2
tol_debias = tolD * (r(:)'*r(:));

% initialize p
p = r;

% main loop
while cont_debias_cg
    
    % calculate A*p =  Rt * R * p ( or if with weights Wt * Rt * R * W * p)
    RWp = AFun(p);
    Ap = AtFun(RWp);
    
    % mask out the zero terms
    Ap(idxZero) = 0;
    
    % calculate alpha for CG
    alpha = rTr / (p(:)'* Ap(:));
    
    % take the step
    x = x + alpha * p;
%     if isNonNegative
%         x(x<0) = 0; % 2017.05.05_Xiang, can't simply do it here. refer to 2013_A Conjugate Gradient Type Method for the Nonnegative Constraints Optimization Problems.: TODO: Apply non negative in a transformed domain
%     end    
    r  = r  - alpha * Ap;
    resid = resid - alpha * RWp; % resid = b - As*xs = resid_prev - As*(xs-xs_prev) = resid_prev - As*(alpha*p) = resid_prev-alpha*RWp
    
    rTrNew = r(:)'*r(:);    
    p = r + rTrNew / rTr * p; % beta = rTrNew / rTr;
    
    rTr = rTrNew;
    
    k = k+1;
    
    objFunVals(k) = 0.5*(resid(:)'*resid(:)) + regWt*regFun(x(:));
    cpuTimes(k) = toc;    
    if ~isempty(xGT)
        [~, SNRs(k)] =  psnr(x, xGT);        
    end    
    % in the debiasing CG phase, always use convergence criterion based on the residual (this is standard for CG)    
    if verbose        
        fprintf('Iteration=%4d, objFunVal=%9.5e, debias resid = %13.8e, convergence = %8.3e, nnz=%7d, nNeg=%7d', k, objFunVals(k),  resid(:)'*resid(:), rTr / tol_debias, compare.nnz(x), compare.nNeg(x))
        if ~isempty(xGT) % xGT
            fprintf(', SNR (x vs xGT) = %4.5edB', SNRs(k));
        end
        fprintf('\n');
    end
    
    cont_debias_cg = (k-debias_start <= minIterD)   ||   ((rTr > tol_debias)  &&  (k-debias_start <= maxIterD));    
end

solver.printAlgorithmEnd('algorithmName', mfilename, 'f_x', objFunVals(end)-regWt*regFun(x), 'regFun_x', regFun(x), 'regWt', regWt, ...
                         'nnz_x', compare.nnz(x), 'CPUTime', cpuTimes(end), 'verbose', verbose);                     

end

