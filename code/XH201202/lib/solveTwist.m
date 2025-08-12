function [x, x_debias, objFunVals, cpuTimes, debias_start, SNRs, sampledHist, regFunVals, maxSVDofA] = solveTwist(varargin)
%[x, x_debias, objFunVals, cpuTimes, debias_start, SNRs, maxSVDofA] = solveTwist(varargin)
%Usage: check solver.parseInput(); for input arguments
%[x, x_debias, objFunVals, cpuTimes, debias_start, SNRs, maxSVDofA] = solveTwist(b, A, regWt, varargin)
%                    [x, x_debias, objFunVals, cpuTimes, debias_start, SNRs, maxSVDofA] = ...
%                         solveTwist(b, AFun, regWt,... 
%                         'AtFun', AtFun, ...                        
%                         'isNonNegative', isNonNegative,...
%                         'denoiseFun', denoiseFun,...
%                         'regFun', regFun,...
%                         'x0', x0,...
%                         'xGT', xGT,...
%                         'maxIterA', maxIterA,...                        
%                         'stopCriterion', stopCriterion,...                                                
%                         'tolA', tolA,...
%                         'isDebias', isDebias, ...
%                         'maxIterD', maxIterD,...
%                         'isMonotone', isMonotone, ...
%                         'verbose', verbose
%                         'lb', lb
%                         'ub', ub
%                          );
% 
% 2020.07: added lower bound and upper bound
%       'lb', -inf, 'ub', +inf as default
%       'isNonNegative', true 
%       is the same, and duplicate as 
%       'lb', 0
%       but we keep 'isNonNegative' which is a dupliate option for backward compatability
%assert((nargin-numel(varargin)) == 3, 'Wrong number of required parameters');

%% Parse Inputs
opts =  solver.parseInput(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=opts.%s;', x, x)), fields(opts)); % exact each field as individual varibles

% alpha, beta, sparse are special reserved names by Matlab, need to explict specify
alpha = opts.alpha;
beta = opts.beta;
sparse = opts.sparse;

% %% Special Case: if regWt is large enough and regFun = L1, the optimal solution is the zero vector
% if solver.isRegWtTooLargeAndL1(opts)
%     warning('The weight of regularizer is too large. Output all zeros!');
%     x = zeros(size(x0));
%     objFunVals(1) = 0.5*(b(:)'*b(:));
%     cpuTime(1) = 0;
%     if ~isempty(xGT)
%         [~, SNRs(1)] = psnr(x, xGT);        
%     end
%     return
% end

%% Initialization
debias_start = 0;
sampledHist = [];
x_debias = [];
SNRs = [];
%%% 1. Common all algorithm initializations: k, x, objFunVals, cpuTimes, SNRs, objFunValGT
% 1.1: k, x(k-2), x(k-1), x(k)
k = 0; %iteration number
xm2 = x0;
xm1 = x0;
x = x0;
% 1.2: objFunVals(k+1), cpuTimes(k+1), SNRs(k+1), objFunValGT
r =  b-AFun(x);
objFunVals(k+1) = objFun(x, r, regWt, regFun);
regFunVals(k+1) = regFun(x); 
tic;  % start the clock
cpuTimes(k+1) = toc;
if ~isempty(xGT)
    [~, SNRs(k+1)] = psnr(x, xGT);
    objFunValGT = objFun(xGT, b-AFun(xGT), regWt, regFun);
else
    objFunValGT = [];
end

%%% 2. Twist only initializations: variables controling first and second order iterations
IST_iters = 0;
TwIST_iters = 0;
rho0 = (1-lambda1/lambdaN)/(1+lambda1/lambdaN);
if isempty(alpha)
    alpha = 2/(1+sqrt(1-rho0^2));
end
if isempty(beta)
    beta  = alpha*2/(lambda1+lambdaN);
end
%--------------------------------------------------------------
%% TwIST iterations: given kth iteration, solve (k+1)th iteration
%--------------------------------------------------------------
solver.printAlgorithmBegin('algorithmName', mfilename, 'verbose', verbose, 'x', x, 'objFunVals', objFunVals, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);    
%while cont_outer
for k = 1:maxIterA    
    % objFunVals(1) is for k=0th initlialize, so objFunVals(k+1) is for kth iteration
    % prepare for next step
    xm2 = xm1;
    xm1 = x;
    
    grad = - AtFun(r); % gradient, we reuse r = b-A*x(k-1) from last step 
    while true
        % IST estimate
        x = denoiseFun(xm1 - grad/maxSVDofA, regWt/maxSVDofA);
%         % Here No need to set x(x<0) = 0 for isNonNegative; % Xiang: 
%         if isNonNegative
%             x(x<0) = 0; % Xiang: TODO: Apply non negative in a transformed domain. Newly added on 2017.04.20
%         end
        if (IST_iters >= 2) || ( TwIST_iters ~= 0)
            % set to zero the past when the present is zero suitable for sparse inducing priors
            if sparse
                mask = (x ~= 0);
                xm1 = xm1.* mask;
                xm2 = xm2.* mask;
            end
            % two-step iteration
            xm2 = (alpha-beta)*xm1 + (1-alpha)*xm2 + beta*x;
            % Lower bound and upper bound.,             
            xm2(xm2 < lb) = lb;
            xm2(xm2 > ub) = ub;            
            % keep isNonNegative for backward compapatbility, though it is may be duplicate of lb
            if isNonNegative
                xm2(xm2 < 0) = 0; % Xiang: TODO: Apply non negative in a transformed domain
            end
            % compute residual
            r = b - AFun(xm2);
            objFunVals(k+1) = objFun(xm2, r, regWt, regFun);
            regFunVals(k+1)  = regFun(x); 
            if (objFunVals(k+1) > objFunVals(k)) && (isMonotone)
                TwIST_iters = 0;  % do a IST iteration if monotonocity fails
            else
                TwIST_iters = TwIST_iters+1; % TwIST iterations
                IST_iters = 0;
                x = xm2;
                if mod(TwIST_iters,10000) == 0
                    maxSVDofA = 0.9*maxSVDofA;
                end
                break;  % break loop while
            end
        else
            % Lower bound and upper bound.,             
            x(x < lb) = lb;
            x(x > ub) = ub;
            % keep isNonNegative for backward compapatbility, though it is may be duplicate of lb
            if isNonNegative
                x(x < 0) = 0; % Xiang: TODO: Apply non negative in a transformed domain
            end
            
            r = b-AFun(x);
            objFunVals(k+1) = objFun(x, r, regWt, regFun);
            if objFunVals(k+1) > objFunVals(k)
                % if monotonicity  fails here  is  because  TODO: set max eig(A'A) to 1 before the solver, for normalization
                % max eig (A'A) > 1. Thus, we increase our guess of max_svs
                maxSVDofA = 2*maxSVDofA;
                if verbose
                    fprintf('Incrementing S=%2.2e\n',maxSVDofA)
                end
                IST_iters = 0;
                TwIST_iters = 0;
            else
                TwIST_iters = TwIST_iters + 1;
                break;  % break loop while
            end
        end
    end
                   
    % store and print out results for current step    
    cpuTimes(k+1) = toc;
    if ~isempty(xGT);  [~, SNRs(k+1)] = psnr(x, xGT); end
    sampledHist = solver.updateSampledHist(sampledHist, x, cpuTimes, k, sampleInterval); % Only keep x history every sampleInterval
    
    criterion = solver.calcConvergeCriterion(stopCriterion, objFunVals(k), objFunVals(k+1), xm1, x, xGT); 
    solver.printAlgorithmIterK('k', k, 'verbose', verbose, 'x', x, 'objFunVals', objFunVals, 'criterionVsTol', criterion/tolA, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);   
    
    if (criterion<=tolA) && (k>minIterA); break;  end    % stop if converged.     
end
solver.printAlgorithmEnd('algorithmName', mfilename, 'cpuTime', cpuTimes(end), 'verbose', verbose, 'xGT', xGT, 'objFunValGT', objFunValGT, 'regWt', regWt, 'regFun', regFun);
%--------------------------------------------------------------
% end of the main loop
%--------------------------------------------------------------

%--------------------------------------------------------------
% we may remove the bias from the l1 penalty, by solve the least-squares problem obtained by omitting the l1 term and fixing the zero elements at zero.
% TODO: Add non negative constraints here
if isDebias
    [x_debias, objFunVals, cpuTimes, debias_start, SNRs] = solveDebias(b, AFun, AtFun, x, xGT, tolD, maxIterD, minIterD, isNonNegative, verbose, k, objFunVals, cpuTimes, SNRs, regWt, regFun, t0);
end
end