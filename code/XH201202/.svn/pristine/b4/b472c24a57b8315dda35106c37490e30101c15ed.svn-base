function [p, gtp, ncg1, dnew, eig_val] = ...
	modlnp_sub (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun)

global free_ind
%---------------------------------------------------------
% this routine performs a preconditioned conjugate-gradient
% iteration to solve the Newton equations for a search
% direction for a truncated-newton algorithm. 
% When the value of the quadratic model is sufficiently 
% reduced, the iteration is terminated.
%---------------------------------------------------------
% parameters
%
% p           - computed search direction
% g           - current gradient
% gv,gz1,v    - scratch vectors
% r           - residual
% d           - diagonal preconditoning matrix
% feval       - value of quadratic function
%------------------------------------------------------------
% initialization
%------------------------------------------------------------
accrcy = argvec(1);
gnorm  = argvec(2);
xnorm  = argvec(3);
%%% --------------------- Reduce the dimension of Hp=-g
bounds=0;
g=g(free_ind);
d=d(free_ind);
%%% -----------------------------------------------------

eig_val = [];

if (maxit == 0);
  p    = -g;
  gtp  = p'*g;
  ncg1 = 1;
  dnew = d;
  if (norm(p)==0);
      disp('MODLNP 01: |p| = 0');
      pause(1);
  end;
  return;
end;
first = 1;
tol   = 1.d-12;  %%%### was 1.d-12  %%%###
qold  = 0;
dnew    = initpc (d, upd1, ireset);
r     = -g;
v     = zeros(size(r));
p     = v;
gtp   = p'*g;

rho  = zeros(maxit+1,1);
beta = zeros(maxit,1);
v_gv = zeros(maxit,1);

rho(1) = norm(r)^2;

%------------------------------------------------------------
% main iteration (conjugate-gradient iterations for Ax = b)
%------------------------------------------------------------
ind = 0;
ncg1 = 0;
for k = 1:maxit
   ncg1 = ncg1 + 1;
   if (bounds); r  = ztime ( r, ipivot); end;
   zk = msolve (r, upd1, ireset, first, d);
   if (bounds); zk = ztime (zk, ipivot); end;
   rz = r'*zk;
   if (rz/gnorm < tol); 
       ind = 80;
       if (norm(p)==0); p = -g; gtp = p'*g; 
       end;
       return; 
   end;
   if (k > 1) 
      beta(k) = rz/rzold;
   else
      beta(k) = 0;
   end;
   v = zk + beta(k)*v;
   if (bounds); v  = ztime ( v, ipivot); end;
   gv = gtims_sub (v, x, g, accrcy, xnorm, sfun);
   if (bounds); gv = ztime (gv, ipivot); end;
   v_gv(k) = v'*gv;
   if (v_gv(k)/gnorm < tol); 
       ind = 50; 
       if (norm(p)==0);
          disp('MODLNP 03: |p| = 0');
          p=-g; gtp=p'*g;
       end;
       return; 
   end;
   ndia3 (dnew, v, gv, r, v_gv(k));
%------------------------------------------------------------
% compute current solution and related vectors
%------------------------------------------------------------
   alpha = rz / v_gv(k);
   p = p + alpha* v;
   r = r - alpha*gv;

   rho(k+1) = norm(r).^2;
   
%------------------------------------------------------------
% Compute the Lanczos spectral estimate.
%------------------------------------------------------------
   eig_val = spectral_estimate (k, rho, v_gv, beta);

%------------------------------------------------------------
% test for convergence
%------------------------------------------------------------
   gtp = p'*g;
   pr = r'*p;
   q = (gtp + pr) / 2;
   qtest = k * (1 - qold/q);
   if (qtest <= 0.5); 
       if (norm(p)==0);
          disp('MODLNP 04: |p| = 0');
          pause(1);
       end;
       return; 
   end;
   qold = q;
%------------------------------------------------------------
% perform cautionary test
%------------------------------------------------------------
   if (gtp > 0); 
       ind = 40; 
       if (norm(p)==0);
          disp('MODLNP 05: |p| = 0');
          pause(1);
       end;
       return; 
   end;
   rzold = rz;
end;
k = k-1;
%------------------------------------------------------------
% terminate algorithm
%------------------------------------------------------------
if (ind == 40)
   p   = p - alpha*v;
end;
%------------------------------------------------------------
if (ind == 50 & k <= 1)
   p = msolve (g, upd1, ireset, first, d);
   p = -p;
   if (bounds); p = ztime (p, ipivot); 
   end;
end;
%------------------------------------------------------------
if (ind == 80 & k <= 1)
   p = -g;
   if (bounds); p = ztime (p, ipivot);
   end;
end;
%------------------------------------------------------------
% store new diagonal preconditioner
%------------------------------------------------------------
gtp  = p'*g;
ncg1 = k + 1;
if (norm(p)==0);
    disp('MODLNP 06: |p| = 0');
    pause(1);
end;
return;

%------------------------------------------------------------------------------
function eig_val = spectral_estimate (k, rho, p_Ap, beta)

%---------------------------------------------
% Estimate the spectrum via Lanczos.

Delta_inv = diag(1./sqrt(rho(1:k)));
B = eye(k) - diag(beta(2:k), 1);
T = Delta_inv*B'*diag(p_Ap(1:k))*B*Delta_inv;
T = 0.5*(T+T');  % This makes sure that T is symmetric!
eig_val = sort(eig(T));
% fprintf('condition number = %e\n', eig_val(end)/eig_val(1));
% fprintf('estimated spectral range: [%e, %e]\n', eig_val(1), eig_val(end));
return;
