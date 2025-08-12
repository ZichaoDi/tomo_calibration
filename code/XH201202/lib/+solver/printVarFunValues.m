function printVarFunValues(x, objFunVal, regWt, regFunVal)
%solver.printVarFunValues(x, objFunVal, regWt, regFunVal): print the variable and function values
%% Print Results
% we use %10.3e instead of %0.3e to handle both positive and negative values with same spaces.
fprintf('Objective function = f(x) + g(x) = %10.3e.\n', objFunVal);
fprintf('    f(x) = 0.5*||Ax-b||_2^2 = %10.3e.\n', objFunVal - regWt*regFunVal); % 
fprintf('    g(x) = regWt*||x||_1 = %.3e * %10.3e = %10.3e.\n', regWt, regFunVal, regWt*regFunVal);
fprintf('Number of non-zero components of x = %d.\n', compare.nnz(x));
end