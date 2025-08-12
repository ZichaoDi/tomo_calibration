function printAlgorithmBegin(varargin)
% solver.printAlgorithmBegin('algorithmName', mfilename, 'verbose', verbose, 'x', x, 'objFunVals', objFunVals, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);
%  show the begining information of a solver

%% Parse Inputs
p = inputParser;  % default CaseSensitive=0, PartialMatching=1
p.addParameter('algorithmName',   '',  @ischar);
p.addParameter('verbose',   true,   @(x) parser.isTrueFalse(x) || isempty(x));
p.addParameter('x',         [],     @(x) isnumeric(x) || isempty(x));
p.addParameter('xGT',       [],     @(x) isnumeric(x) || isempty(x));
p.addParameter('objFunVals',[],     @(x) parser.isNumericScalarOrVector(x) || isempty(x));  %
p.addParameter('objFunValGT', [],   @(x) parser.isNumericScalar(x) || isempty(x));
p.addParameter('regWt',       [],     @(x) parser.isPositiveScalar(x) || isempty(x));
p.addParameter('regFun',    [],     @(x) isa(x, 'function_handle')); 
p.parse(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=p.Results.%s;', x, x)), fields(p.Results)); % extract variables from structure, or just use in.variable after set in =p.Results;

%% Print Results
if verbose            
    % 1. Print starting mark
	fprintf('************************************************************************************************\n');    
    fprintf('****************************     Start algorithm ''%s'':     ****************************\n', algorithmName); 
    % 1. Print definition of f,g,regWt
    fprintf('**** Objective function = f(x) + g(x) = 0.5*||Ax-b||_2^2 + regWt*regFun(x), where regWt = %.3e. ****\n', regWt);
    % 3. Print initial step
    solver.printAlgorithmIterK('k', 0, 'verbose', verbose, 'x', x, 'objFunVals', objFunVals, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);
    %solver.printVarFunValues(x, objFunVal, regWt, regFun(x))
end

end

