 function [ipivot, ierror, xnew] = crash (x, low, up)
%---------------------------------------------------------
% this initializes the constraint information, and
% ensures that the initial point satisfies 
%      low <= x <= up.
% the constraints are checked for consistency.
%---------------------------------------------------------
ierror = 0;
if (any(low > up)); ierror = - max(find(low > up)); end;
x      = max (low, x);
x      = min ( up, x);
xnew   = x;
ind    = find(abs(low-x)<10*eps);
ipivot = zeros(size(x));
if (length(ind) > 0); ipivot(ind) =  -ones(size(x(ind))); end;
ind    = find(abs(up-x)<10*eps);
if (length(ind) > 0); ipivot(ind) =   ones(size(x(ind))); end;
ind    = find(low == up);
if (length(ind) > 0); ipivot(ind) = 2*ones(size(x(ind))); end;
