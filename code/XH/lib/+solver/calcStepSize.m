function alpha = calcStepSize(x, dx, maxStepSize)
%solver.calcStepSize()  caculate the largest stepsize alpha in (0, 1] that makes x + alpha*dx >= (1- maxStepSize)*x
% maxStepsize is set to 1

% use isAllowFullNewtonStep = true for XH paper, false for 2014_Fountoulakis, Gondzio_Matrix-free Interior Point Method for Compressed Sensing Problems
% There is no reason to NOT allow full newton step as in Gondzio as long as the (x,s)>0 and not too small
isAllowFullNewtonStep = true; % true or false; 
if isAllowFullNewtonStep
    maxStepSize = 1; 
else
    if ~exist('maxStepSize', 'var') || isempty(maxStepSize)
        maxStepSize = 0.999; % was using 0.9999 before 2018.02.07,  0.99~0.9999 is ok; % 0.999999 is too big, 0.9 is too small,
    end
    assert(maxStepSize<=1 && maxStepSize>0, sprintf('maxStepSize need to bewtween (0, 1], but it is %.8e', maxStepSize));
end

gamma = 0.999;

iNeg = (dx < 0);
if nnz(iNeg) > 0
    alphaMax = min(x(iNeg)./ (-dx(iNeg))); % alpha = min(1, min(x(iNeg)./ (-dx(iNeg)))); % add 0.999 to make non negtaive
else
    alphaMax = Inf;    % all dx >= 0, any positive step size makes x + alpha*dx > 0
end

alpha = min(maxStepSize, gamma*alphaMax);

end
