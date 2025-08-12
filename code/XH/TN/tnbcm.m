function [xstar, f, g, ierror, varargout] = tnbcm (x, sfun, low, up, maxiter)
%---------------------------------------------------------
% this routine solves the optimization problem
%
%   minimize     f(x)
%   subject to   low <= x <= up
%
% parameters:
%
% ierror  <-  error code
%             ( 0 => normal return
%             ( 2 => more than maxfun evaluations
%             ( 3 => line search failed (may not be serious)
%             (-1 => error in input parameters
% x        -> initial estimate of the solution; 
% sfun     -> function routine: [f,g] = sfun(x)
% xstar   <-  the computed solution.
% g       <-  final value of the gradient
% f       <-  final value of the objective function
% low, up  -> lower and upper bounds on the variables
% maxiter  -> upper limit on number of nonlinear iterations
%
% this routine sets up all the parameters for lmqnbc:
%
% maxfun - maximum allowable number of function evaluations
% stepmx - maximum allowable step in the linesearch
% accrcy - accuracy of computed function values
% maxit  - maximum number of inner iterations per step
%---------------------------------------------------------
global tn_max_maxit tn_maxfun
%---------------------------------------------------------
if (isempty(tn_max_maxit))
  max_maxit=50;
else
  max_maxit = tn_max_maxit;
end;
n      = length(x);
maxit  = 1 + round((n+1)/2);
if (exist('max_maxit'))
  maxit  = min(max_maxit, maxit);
else
  maxit  = min(50, maxit);
end;
if (isempty(tn_maxfun))
  maxfun = 150*n;
else
  maxfun = tn_maxfun;
end
stepmx = 10;
accrcy = 100*eps;
%---------------------------------------------------------
[xstar, f, g, ierror, eig_val] = ...
   lmqnbcm1 (x, sfun, low, up, maxiter, maxit, maxfun, stepmx, accrcy);

if (nargout == 5)
  varargout{1} = eig_val;
end
