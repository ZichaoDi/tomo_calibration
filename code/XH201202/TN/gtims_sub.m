function gv = gtims_sub (v, x, g, accrcy, xnorm, sfun)
global free_ind
%%%%%%%%%%%%%%%%%##################
%---------------------------------------------------------
% compute the product of the Hessian times the vector v;
% store result in the vector gv 
% (finite-difference version)
%---------------------------------------------------------
%%%%%%%%%%%%%%%%%##################
%%##  delta   = sqrt(accrcy)*(1 + xnorm);
delta   = sqrt(accrcy)*(1 + xnorm)/norm(v);
v_full=zeros(size(x));v_full(free_ind)=v;
hg      = x + delta*v_full;
%%%%%%%%%%%%%%%%%##################
[~, gv] = feval (sfun, hg);
gv      = (gv(free_ind) - g)/delta;
