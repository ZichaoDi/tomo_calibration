function opts = parseInput(varargin)
%solver.parseInput: parse input parameters for a optimization solver
% Input Cases:
% (0) All inputs are name/valued pairs: varargin = {'name1', value1, 'name2', value2, ...}.
% (1) Inputs is a single structure, varargin = {opts}. it is already parsed options, just return opts = varargin.
% (2) First input is mfm object, rest are name/valued pairs: varargin = {mfm, 'name1', value1, 'name2', value2, ...}.
%    NOTE: (1) & (2) are distinguishable because (1) is a single structure, (2) is a single object followed by name/valued pairs
% (3) First three inputs are b, A, rest are name/valued pairs: varargin = {b, A, 'name1', value1, 'name2', value2, ...}.
% see SOLVER.opts for each meaning of name/value pair

%% (1)/(2)/(3) In those cases, the first input is not string
if ~isstr(varargin{1})  % Matlab suggests to use ischar(s) instead of isstr(s) for check if it s is a string. No, I like isstr for string
    if isstruct(varargin{1})
        % (1) Inputs is a single structure object, varargin = {opts}.
        assert(nargin == 1, 'Unkown input formats! When first input is a structure, it should not be followed by more inputs!');
        opts = varargin{1};
    elseif isobject(varargin{1})
        % (2) varargin = {mfm, 'name1', value1, 'name2', value2, ...}.
        % b/AFun/AtFun/xGT/diagATA/Vk can be obtaind from mfm class
        mfm = varargin{1};
        opts = solver.parseInput('b', mfm.g, 'AFun', mfm.AFun, 'AtFun', mfm.AtFun, 'xGT', mfm.fGT, 'diagATA', mfm.diagATA, 'Vk', mfm.Vk, varargin{2:end});
    elseif isnumeric(varargin{1})
        % (3) varargin = {b, A, 'name1', value1, 'name2', value2, ...}.
        assert(isnumeric(varargin{2}), 'Unkown input formats! When first input is numeric (b), second input should also be numeric (A)');
        opts = solver.parseInput('b', varargin{1}, 'AFun', varargin{2}, varargin{3:end});
    else
        error('Unkown input formats!');
    end
    return
end


%% (3) varargin = {}
% (0) All inputs are name/valued pairs: varargin = {'name1', value1, 'name2', value2, ...}.

%% (0) varargin = {'name1', value1, 'name2', value2, ...}.
%% L1 norm and soft threshold function
L1Fun = @(x) sum(abs(x(:)));                % L1 norm,
softFun = @(x,T) sign(x).*max(abs(x)-T, 0); % soft threshold for real. If complex, b = max(abs(x) - T, 0); b = b./(b+T) .* x;

%% Parse
p = inputParser;  % default CaseSensitive=0, PartialMatching=1
isTrueFalseOrEmpty = @(x) isempty(x) || isscalar(x) && ismember(x, [true, false, 1, 0]); % true/false/1/0
isPositiveFloatOrEmpty     = @(x) isempty(x) || isnumeric(x) && isscalar(x) && (x>0);
isPositiveIntegerOrEmpty   = @(x) isempty(x) || isnumeric(x) && isscalar(x) && (x>0) && (x==round(x));

% New SOlVE class only parameters
%p.addParameter('solType', 'twist',  @(x) isempty(x) || ischar(x) && ismember(lower(x), {'twist', 'fista', 'gpsr', 'spiral', 'twist_nonnegative'}));  % choose from {'twist', 'fista', 'spiral'}
p.addParameter('solType', 'twist');

p.addParameter('b',     [],     @(x) isnumeric(x));
p.addParameter('bSz',   []);  % Inferred from size(b)
p.addParameter('AFun',    [],     @(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x)); % could be either a matrix or a function handle, or []. If it is a matrix, pass it to AMatrix, and reassign it as handle
p.addParameter('AtFun',   [],     @(x) isa(x, 'function_handle') || isempty(x)); % must be function handle or []
p.addParameter('xGT',	[],     @(x) isempty(x) || isnumeric(x));  % numeric array
p.addParameter('x0',	'Atb',  @(x) isempty(x) || isnumeric(x) || ischar(x) && ismember(lower(x), {'zero', 'zeros', 'random', 'rand', lower('Atb')}));  % numeric array or string choice
p.addParameter('AMatrix',     [],     @(x) isempty(x) || isnumeric(x)); % matrix form or empty (not specified)
p.addParameter('xSz',   []);  % Inferred from AtFun(b), or speficied when xSz can't be inferred: if AFun is matrix A and both x0/xGT are not numeric array.

p.addParameter('regWt',         1,      isPositiveFloatOrEmpty);
p.addParameter('regFun',        L1Fun,  @(x) isempty(x) || isa(x, 'function_handle') || ischar(x) && ismember(lower(x), {'l0', 'l1', 'tv', 'wavelet'}));  % default L1 norm
p.addParameter('objFun',        [],     @(x) isempty(x) || isa(x, 'function_handle'));
p.addParameter('denoiseFun',    softFun, @(x) isempty(x) || isa(x, 'function_handle'));  % default soft threshold
p.addParameter('Transform',     false,  isTrueFalseOrEmpty); % TODO, wavelet, DCT, TV etc

% Display Parameters
p.addParameter('verbose',   true,  isTrueFalseOrEmpty);


% General Parameters: Used for both TwIST and FISTA
%            p.addParameter('regFun',   L1Fun,      @(x) isempty(x) || isa(x, 'function_handle'));  % default L1 norm
%            p.addParameter('denoiseFun',    softFun,    @(x) isempty(x) || isa(x, 'function_handle'));  % default soft threshold
%            p.addParameter('AtFun',   [],         @(x) isempty(x) || isa(x, 'function_handle'));  % default soft threshold
p.addParameter('stopCriterion',	1,  @(x) isempty(x) || isstr(x) || parser.isNumericScalar(x));
p.addParameter('tolA',          1e-6,   isPositiveFloatOrEmpty); % or 1e-2
p.addParameter('maxIterA',      100,    isPositiveIntegerOrEmpty);
p.addParameter('minIterA',      5,      isPositiveIntegerOrEmpty);
p.addParameter('isNonNegative',	true,  isTrueFalseOrEmpty); % TODO, NEW parameter whether x is non-negative or not
% Debais Parameters (optional step)
p.addParameter('isDebias',  false,  isTrueFalseOrEmpty);
p.addParameter('tolD',      1e-7,   isPositiveFloatOrEmpty); % or 1e-3
p.addParameter('maxIterD',  200,     isPositiveIntegerOrEmpty); % or 10
p.addParameter('minIterD',  5,     isPositiveIntegerOrEmpty);
p.addParameter('sampleInterval',  3600,     @(x) isempty(x) || isnumeric(x) && all(x>0) || iscell(x)); % 1 hour = 3600 seconds.

% Parameters for precondtioner
p.addParameter('preconditioner',  [],   @(x) isempty(x) || isa(x, 'function_handle') || isnumeric(x) || isstr(x) || iscell(x)&&numel(x)==2); % []/matrix M/ {}
p.addParameter('diagATA',   [],   @(x) isempty(x) || isnumeric(x));
p.addParameter('Vk',        [],   @(x) isempty(x) || isnumeric(x));

% IP Only Parameters
p.addParameter('muCG',          []);

% TwIST Only Parameters
p.addParameter('isMonotone',	true,	isTrueFalseOrEmpty);
p.addParameter('lambda1',       1e-4,   isPositiveFloatOrEmpty);
p.addParameter('lambdaN',       1,      isPositiveFloatOrEmpty); % or 1e-2
p.addParameter('sparse',        true,   isTrueFalseOrEmpty);
p.addParameter('alpha',         [],   isPositiveFloatOrEmpty);
p.addParameter('beta',          [],   isPositiveFloatOrEmpty);
p.addParameter('maxSVDofA',     1e-8,   @(x) isnumeric(x)); % 2018.11.12 added: max singular value of A, default 1e-8 (in original TwIST paper is 1)


% FISTA Only Parameters
p.addParameter('Lipschitz0',    1,     isPositiveFloatOrEmpty);  % a value smaller than Lipschitz constant

% TODO: put the normalizer in the solver, not in the model!!!
p.parse(varargin{:});
opts = p.Results;

%% Further Processing
% Whenever we provide the full matrix via opts.AMatrix or opts.AFun, we assign/reassign AFun/AtFun as function handle using A. Then no FFT fast operation of A*x and A'*b.
if (~isempty(opts.AMatrix) && isnumeric(opts.AMatrix)) || (~isempty(opts.AFun) && isnumeric(opts.AFun))
    if isempty(opts.AMatrix)
        opts.AMatrix = opts.AFun;
    end
    opts.AFun  = @(x) reshape(opts.AMatrix * x(:), size(opts.b));
    opts.AtFun = @(b) reshape(opts.AMatrix' * b(:), opts.xSz);
end

% parse if opts.regFun is a string specifies function choice
if ~isempty(opts.regFun) && ischar(opts.regFun)
    switch lower(opts.regFun)
        case lower('L0')
            opts.regFun = @(f) l0norm(f);    % L0 norm: could be bad and not converging
            opts.denoiseFun = @(f,th) hard(f,th); % hard denoise may only good for L0-norm, sometimes is too sensitive and hard to get good result
        case lower('L1')
            opts.regFun = @(f) sum(abs(f(:)));  % L1 norm: cube looked like sphere for 20160421_1134_simulated_sigma=0.00_nVoxel=98598_L1NormOriginal
            %opts.denoiseFun = @(f,th) soft(f,th); %opts.denoiseFun = @(f,t) sign(f).*max(abs(f)-th,0);
            opts.denoiseFun = @(f,th) sign(f).*max(abs(f)-th,0);
        case lower('TV')
            % Total variation, could be 3D or 2D
            tv_iters = 5;            
            opts.regFun = @(f)TVnorm3D(f, [], 'isotropic'); %opts.regFun = @(f)TVnorm3D(f, mfm.voxSz, 'isotropic'); %istropic or anistropic
            opts.denoiseFun = @(f,th)TVdenoise3D(f, 2/th, tv_iters); % TODO: add voxSz there?
        case lower('wavelet')
            optw.wname = 'db4'; optw.wlev = 3; % 'db4', 'sym8','coif4'; wavelet name
            opts.regFun = @(f) L1Wavelet3D(f,optw); % Wavelet: not much difference compare with TVnorm3D
    end
end


opts.bSz = size(opts.b);
Atb = opts.AtFun(opts.b);
if isempty(opts.xSz)
    opts.xSz = size(Atb);
else
    assert(all(opts.xSz == size(Atb)), 'Inconsistent size of x!');
end
opts.x0 = solver.calcInitalX0(opts.x0, opts.AFun, opts.AtFun, opts.b, Atb);
if ~isempty(opts.x0)
    assert(isscalar(opts.regWt) || all(size(opts.regWt)==opts.xSz), 'Parameter regWt has wrong dimensions; it should be scalar or size(x)');
    assert(isempty(opts.xGT) || all(size(opts.xGT)==opts.xSz), 'xGT has incompatible size: need to be [] or compatiable with A');
    assert(isempty(opts.x0) || all(size(opts.x0)==opts.xSz), 'x0 has incompatible size: need to be [] or compatiable with A');
end

if isempty(opts.objFun)
    opts.objFun = @(x, resid, regWt, regFun) 0.5*(resid(:)'*resid(:)) + regWt*regFun(x);
end

if exist('A_operator.m', 'file') % check if the A_operator class is defined, don't use exist('@A_operator', 'dir') as even we put @A_operator into a folder in search path, matlab may claim @A_operator not exist if we search from a different folder, while A_operator.m always could be found.
    opts.A = A_operator(@(x) reshape(opts.AFun(reshape(x, opts.xSz)),  [], 1),...
                        @(b) reshape(opts.AtFun(reshape(b, opts.bSz)), [], 1));
end


end
