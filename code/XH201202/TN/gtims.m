function gv = gtims (v, x, g, accrcy, xnorm, sfun)

%%%%%%%%%%%%%%%%%##################
%---------------------------------------------------------
% compute the product of the Hessian times the vector v;
% store result in the vector gv 
% (finite-difference version)
%---------------------------------------------------------
%%%%%%%%%%%%%%%%%##################
%%##  delta   = sqrt(accrcy)*(1 + xnorm);
delta   = sqrt(accrcy)*(1 + xnorm)/norm(v);
hg      = x + delta*v;
%%%%%%%%%%%%%%%%%##################
[~, gv] = feval (sfun, hg);
gv      = (gv - g)/delta;
