function x1 = ztime(x,ipivot)
%---------------------------------------------------------
% this routine multiplies the vector x by the 
% constraint matrix z
%---------------------------------------------------------
ind    = find(ipivot~=0);
x(ind) = zeros(size(ind));
x1     = x;
