function printAlgorithmIterK(varargin)
%solver.printAlgorithmIterK('k', k, 'verbose', verbose, 'x', x, 'objFunVals', objFunVals, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);
%  show the information after kth iteration of a solver
% Since k start with 0 for the initial iteration, so objFunVals(k+1) etc are for the kth iteration
% k = 0:            Initial
% k = 1, 2, ..:     kth iteration
% k = []/inf/nan:   Grouth truth
% 2017.09.17
% e.g.
% Iter# =   1: objFunVal = 5.69689627e-02, criterion/tol = 6.8e+14, f = 5.673e-02, g = 2.407e-04, |x|= 2.407e+00, nnz =539014, nNeg =    0, psnrSparse = 1.25e+01dB

%% Parse Inputs
p = inputParser;  % default CaseSensitive=0, PartialMatching=1
p.addParameter('k',         [],     @(x) parser.isPositiveInteger(x+1) || isempty(x) || isinf(x) || isnan(x));
p.addParameter('verbose',   true,   @(x) parser.isTrueFalse(x) || isempty(x));
p.addParameter('x',         [],     @(x) isnumeric(x) || isempty(x));
p.addParameter('xGT',       [],     @(x) isnumeric(x) || isempty(x));
p.addParameter('objFunVals',[],     @(x) parser.isNumericScalarOrVector(x) || isempty(x));  
p.addParameter('objFunValGT', [],   @(x) parser.isNumericScalar(x) || isempty(x));
p.addParameter('regWt',       [],     @(x) parser.isPositiveScalar(x) || isempty(x));
%p.addParameter('criterionVsTol',       [],     @(x) parser.isNonNegativeScalar(x) || isempty(x) || isnan(criterionVsTol)); % when its very small Nan
p.addParameter('criterionVsTol',       []); % when its very small Nan
p.addParameter('regFun',    [],     @(x) isa(x, 'function_handle')); 
% 2017.07.30_added duality
p.addParameter('dualFunVals',[],     @(x) parser.isNumericScalarOrVector(x) || isempty(x));  
p.parse(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=p.Results.%s;', x, x)), fields(p.Results)); % extract variables from structure, or just use in.variable after set in =p.Results;

%% Print Results
% we use %12.5e instead of %0.3e to handle both positive and negative values with same spaces.
if verbose    
    
    %% Set iterNoStr, objFunValChangeStr, 
    getIterNoStr = @(x) sprintf('Iter# =%4d', x);  % show iteration number only for k >0
    iterNoStrInitial = 'Initialize';   % k = 0, 
    iterNoStrGT      = 'GroundTruth';   % k = []/inf/nan
    iterNoStrLen = max([numel(getIterNoStr(1)), numel('Initialize'), numel('GroundTruth')]);
    
    
    %% Converge Criterion
    criterionVsTolStrFun = @(x) sprintf(', criterion/tol =%8.1e', x);     % show decrease only for k > 0    
    criterionVsTolStrSpaces = repmat(' ', 1, numel(criterionVsTolStrFun(0.1)));  % show spaces for other cases k = 0/[]/inf/nan
    
    if isempty(k) || isinf(k) || isnan(k)
        % Print for Ground Truth        
        curObjFunVal = objFunValGT;        
        iterNoStr = iterNoStrGT;
        criterionVsTolStr = criterionVsTolStrSpaces;
    else
        % Print for Iteration k
        curObjFunVal = objFunVals(k+1);
        if k==0
            % Intialize k = 0            
            iterNoStr = iterNoStrInitial;
            criterionVsTolStr = criterionVsTolStrSpaces;
        else
            % Iteration k>0
            iterNoStr = getIterNoStr(k);            
            preObjFunVal = objFunVals(k);
            criterionVsTolStr = criterionVsTolStrFun(criterionVsTol);
        end        
    end
    % left-justify the iterNoStr with fixed width of iterNoStrLen.  '-' means left-justify the input, '*' sets the filed width
    iterNoStr = sprintf('%-*s', iterNoStrLen, iterNoStr); 
    
    regFun_x = regFun(x);
    gx = regWt*regFun_x;
    fx = curObjFunVal - gx;
    
    %% Set recErrorStr with either psnrSparse (if the function to compute psnrSparse avaliable) or regular snr
    if ~isempty(xGT) 
        if exist('difference', 'file')
            [~, psnrSparse] = difference(x, xGT);
            recErrorStr = sprintf(', psnrSparse =%9.4edB', psnrSparse);
        else
            [~, snr_] = psnr(x, xGT);        
            recErrorStr = sprintf(', SNR =%9.4edB', snr_);
        end
    else
        %recErrorStr = repmat(' ', 1, numel(sprintf(', SNR =%9.2edB', 1.0)));
        recErrorStr = '';
    end
    
    %% Set duality values 2017.07.30_added
    dualStr = '';
    if ~isempty(dualFunVals)
        dualStr = sprintf(', dualFunVal =%15.8e, dualGap =%8.1e', dualFunVals(k+1), objFunVals(k+1)-dualFunVals(k+1));
    end
    
    %% Print
    fprintf('%s: objFunVal =%15.8e%s%s, f =%10.3e, g =%10.3e, |x|=%10.3e, nnz =%5d (threshold 1e-3), nNeg =%5d%s\n',...
            iterNoStr, curObjFunVal, dualStr, criterionVsTolStr, fx, gx, regFun_x, compare.nnz(x, 1e-3), compare.nNeg(x), recErrorStr);
    %fprintf('Iteration=%4d, objFunVals=%9.5e, criterion=%7.3e, nnz=%7d, nNeg=%7d', k, objFunVals(k), criterion/tolA, compare.nnz(x), compare.nNeg(x))    
end
end
