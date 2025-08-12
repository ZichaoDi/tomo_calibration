function printAlgorithmEnd(varargin)
%solver.printAlgorithmEnd('algorithmName', mfilename, 'cpuTime', cpuTimes(end), 'verbose', verbose, 'xGT', xGT, 'objFunValGT', objFunValGT, 'regWt', regWt, 'regFun', regFun);
%   show the final result of a solver

%% Parse Inputs
p = inputParser;  % default CaseSensitive=0, PartialMatching=1
p.addParameter('algorithmName',   '',  @ischar);
p.addParameter('xGT',       [],     @(x) isnumeric(x) || isempty(x));
p.addParameter('objFunValGT', [],   @(x) parser.isNumericScalar(x) || isempty(x));
p.addParameter('regWt',       [],     @(x) parser.isPositiveScalar(x) || isempty(x));
p.addParameter('regFun',    [],     @(x) isa(x, 'function_handle')); 
p.addParameter('cpuTime',   [],     @(x) parser.isNonNegativeScalar(x) || isempty(x));
p.addParameter('verbose',   true,   @(x) parser.isTrueFalse(x) || isempty(x));
p.parse(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=p.Results.%s;', x, x)), fields(p.Results)); % extract variables from structure, or just use in.variable after set in =p.Results;

%% Print Results
if verbose    
    % 1. Print Grouth Truth when xGT exists
    if ~isempty(xGT)
        solver.printAlgorithmIterK('k', [], 'verbose', verbose, 'x', xGT, 'objFunValGT', objFunValGT, 'regWt', regWt, 'regFun', regFun, 'xGT', xGT);
    end
    % 2. Print total CPU time
    fprintf('Total CPU time = %.3e seconds.\n', cpuTime);
    % 3. Print finishing mark
    fprintf('****************************     Finish algorithm ''%s'':     ****************************\n', algorithmName);
    fprintf('************************************************************************************************\n');
end

end

%%   Old code with redundent information
%     solver.printVarFunValues(x, objFunVal, regWt, regFun(x))
%     if ~isempty(xGT)
%         fprintf('---------------- Compare with Ground Truth Below: ---------------- \n');
%         solver.printVarFunValues(xGT, objFunValGT, regWt, regFun(xGT))
%     end