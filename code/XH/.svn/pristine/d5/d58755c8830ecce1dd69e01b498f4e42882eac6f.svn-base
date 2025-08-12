function z = ssbfgs (s, hv, hy, ys, yhy, vs, vhy)
%---------------------------------------------------------
% self-scaled bfgs quasi-Newton update 
% (used by preconditioner)
%---------------------------------------------------------
delta = (1 + yhy/ys)*vs/ys - vhy/ys;
beta  = -vs/ys;
z     = hv + delta*s + beta*hy;
