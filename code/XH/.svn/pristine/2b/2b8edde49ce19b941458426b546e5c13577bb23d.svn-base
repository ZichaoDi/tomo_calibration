function e = ndia3 (e, v, gv, r, vgv)
%---------------------------------------------------------
% update the preconditioning matrix based on a diagonal 
% version of the bfgs quasi-newton update.
%---------------------------------------------------------
tol    = 1.d-6;
vr     = v'*r;
if (abs(vr)>tol & abs(vgv)>tol);
  e      = e  - (r.*r)/vr + (gv.*gv)/vgv;
  ind    = find(e < tol);
  if (length(ind) > 0); e(ind) = ones(size(e(ind))); end;
end;
