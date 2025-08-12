function  criterion = calcConvergeCriterion(stopCriterion, objFunValPre, objFunValCur, xPre, xCur, xGT)
%, , iter,miniter,stopCriterion,tolerance, dx, x, cpuTimes, objFunVals)
% =====================================
% = Termination Criteria Computation: =
% =====================================
% outside control loop is that for k = 1:maxIter;  if (k<minIter)
%
%
% another possiblility is to create a function
% function [isConverged, criterion] = checkConvergence(stopCriterion, tolA, objFunValPre, objFunValCur, xPre, xCur, xGT)
% where isConverged  = criterion < tolA
% edge conditions, where we 
% if iter < miniter
%     isStop = false;
%     return
% end
% 
% if cpuTime > maxCpuTime
%     isStop = true;
%     status = 'maximum CPU time reached';
%     return        
% end

switch lower(stopCriterion)
    case {1, 'objective', lower('objectiveChange'), lower('objChange')} 
        % relative changes of the objective function value (assume not 0): typically decrease, but increase also used in some cases).
        criterion =  abs(objFunValPre - objFunValCur) / abs(objFunValPre); %abs(curObjFunVal-objFunVals(k+1))/objFunVals(k+1); % objFunVals(k+1) is the prev_objective;
    case {2,  'x', 'variable', lower('xChange')}
        % relative changes of the variable x (estimate)
        criterion = norm(xCur(:) - xPre(:)) / norm(xPre(:));
    case {3, 'sparsity', lower('nonzeroChange')}
        % absolute change of the sparsity (number of non-zero components) of the variable x
        [~, nzMaskPre] = compare.nnz(xPre);
        [~, nzMaskCur] = compare.nnz(xCur);        
        criterion = sum(nzMaskCur(:) ~= nzMaskPre(:));        
    case {4, lower('GT'), lower('GroundTruth'), lower('differenceFromXGT')}
        assert(exist('xGT', 'var') && ~isempty(xGT), 'xGT should be exit non-empty values!');
        criterion = norm(x(:)-xGT(:)) / norm(xGT(:));
    case {5, 'derivative', lower('smallDerivative')}
        % derivative value
        error('dervivative magnitude condition not implemented yet!');
    case {6, 'dual', lower('dualGap'), lower('dualityGap')}
        % complementarity condition
        error('Duality gap omplementarity condition not implemented yet!');
        %criterion = abs(sum(abs(xCur(:)))- lambda.'*b);
        % Norm of lagrangian gradient        
    otherwise
        error('Unknown stopping criterion!');
end
end