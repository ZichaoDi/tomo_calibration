function [xstar, f, g, ierror, eig_val] = ...
    lmqnbcm (x, sfun, low, up, maxiter, maxit, maxfun, stepmx, accrcy)
%---------------------------------------------------------
% This is a bounds-constrained truncated-newton method.
% The truncated-newton method is preconditioned by a
% limited-memory quasi-newton method (computed by
% this routine) with a diagonal scaling (routine ndia3).
% For further details, see routine tnbc.
%---------------------------------------------------------
global sk yk sr yr yksk yrsr
global NF N current_n lambdaind 
global err0 WS
%---------------------------------------------------------
% check that initial x is feasible and that the bounds
% are consistent
%---------------------------------------------------------
[ipivot, ierror, x] = crash(x, low, up);
if (ierror ~= 0);
    fprintf(1,'LMQNBC: terminating (no feasible point)');
    f = 0;
    g = zeros(size(x));
    xstar = x;
    exit;
end;

%---------------------------------------------------------
% initialize variables, parameters, and constants
%---------------------------------------------------------
fprintf(1,'  it     nf     cg           f           |g|\n')
upd1   = 1;
ncg    = 0;
conv   = 0;
xnorm  = norm(x,'inf');
ierror = 0;
nind   = find(N==current_n);
eig_val{1} = [];
if (stepmx < sqrt(accrcy) | maxfun < 1);
    ierror = -1;
    xstar = x;
    return;
end;
%---------------------------------------------------------
% compute initial function value and related information
%---------------------------------------------------------
[f, g] = feval (sfun, x);
oldf   = f;
gnorm  = norm(g,'inf');
nf     = 1;
nit    = 0;
flast  = f;
%---------------------------------------------------------
% Test if Lagrange multipliers are non-negative.
% Because the constraints are only bounds, the Lagrange
% multipliers are components of the gradient.
% Then form the projected gradient.
%---------------------------------------------------------
ind = find((ipivot ~= 2) & (ipivot.*g>0));
if (~isempty(ind));
    ipivot(ind) = zeros(length(ind),1);
end;

%%%%%%%%%%%%%#######################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%find multiplier index
lamind=find(ipivot ~= 0 & ipivot ~= 2);
toll = 10 * eps * (abs(low) + ones(size(x)));
subl=find(x-low<toll);
lambdaind1=intersect(lamind,subl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tolu = 10 * eps * (abs(up) + ones(size(x)));
subu=find(up-x<tolu);
lambdaindu=intersect(lamind,subu);
lambdaind=struct('low',lambdaind1,'up',lambdaindu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = ztime (g, ipivot);
gnorm = norm(g,'inf');
fprintf(1,'%4i   %4i   %4i   % .8e   %.1e\n', ...
    nit, nf, ncg, f, gnorm)
%---------------------------------------------------------
% check if the initial point is a local minimum.
%---------------------------------------------------------
ftest = 1 + abs(f);
if (gnorm < .01*sqrt(eps)*ftest);
    xstar = x;
    NF(1,nind) = NF(1,nind) + nit;
    NF(2,nind) = NF(2,nind) + nf;
    NF(3,nind) = NF(3,nind) + ncg;
    return;
end;
%---------------------------------------------------------
% set initial values to other parameters
%---------------------------------------------------------
n      = length(x);
icycle = n-1;
ireset = 0;
bounds = 1;
difnew = 0;
epsred = .05;
fkeep  = f;
d      = ones(n,1);
%---------------------------------------------------------
% ..........main iterative loop..........
%---------------------------------------------------------
% compute the search direction
%---------------------------------------------------------
argvec = [accrcy gnorm xnorm];
[p, gtp, ncg1, d, eig_val{nit+1}] = ...
    modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
%%##################################################################
tembl=find(abs(x-low)<eps*10);
if (~isempty(tembl))
    for i=1:length(tembl)
        if (p(tembl(i))<0)
            ipivot(tembl(i))=-1;
            p(tembl(i))=0;
        end
    end
end
tembu=find(abs(x-up)<eps*10);
if (~isempty(tembu))
    for i=1:length(tembu)
        if (p(tembu(i))>0)
            ipivot(tembu(i))=1;
            p(tembu(i))=0;
        end
    end
end
%##################################################################

ncg = ncg + ncg1;
while (~conv);
    oldg = g;
    pnorm = norm(p, 'inf');
    oldf = f;
    %---------------------------------------------------------
    % line search
    %---------------------------------------------------------
    pe = pnorm + eps;
    spe = stpmax (stepmx, pe, x, p, ipivot, low, up);
    alpha = step1 (f, gtp, spe);
    alpha0 = alpha;
    PieceLinear=1;
    if(PieceLinear)
        [x, f, g, nf1, ierror, alpha] = lin_proj (p, x, f, g, alpha0, sfun, low, up);
    else
        [x, f, g, nf1, ierror, alpha] = lin1 (p, x, f, alpha0, g, sfun);
    end
    if(current_n==N(1))
    ErrIter=norm(x-WS(:))/err0
    end
    %---------------------------------------------------------
    %## MODIFIED: ADDED THIRD CONDITION TO 'IF' STATEMENT ##
    if (alpha == 0 && alpha0 ~= 0);
        format long e
        fprintf('Error in Line Search (lmqnbcm.m)\n');
        fprintf('    ncg1   = %d\n',ncg1);
        fprintf('    alpha  = %12.8f\n', alpha);
        fprintf('    alpha0 = %12.8f\n', alpha0);
        fprintf('    g''p    = %12.4e\n', gtp);
        fprintf('    |g|     = %12.4e\n', norm(g));
        fprintf('    |p|     = %12.4e\n', norm(p));
        ierror = 2;
        xstar = x;
        if isempty(ncg); ncg=0; end;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    nf  = nf  + nf1;
    nit = nit +   1;
    %---------------------------------------------------------
    % update active set, if appropriate
    %---------------------------------------------------------
    newcon = 0;
    if (abs(alpha-spe) <= 10*eps);
        ipivot0 = ipivot;
        newcon = 1;
        [ipivot, flast] = modz (x, p, ipivot, low, up, flast, f, alpha);
        dipiv = norm(ipivot-ipivot0);
        if (dipiv > 0); ierror = 0; end;
    end;
    if (ierror == 3);
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    %---------------------------------------------------------
    % stop if more than maxfun evaluations have been made
    %---------------------------------------------------------
    if (nf > maxfun);
        ierror = 2;
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    %---------------------------------------------------------
    % set up for convergence and resetting tests
    %---------------------------------------------------------
    difold = difnew;
    difnew = oldf - f;
    if (icycle == 1);
        if (difnew >  2*difold); epsred =  2*epsred; end;
        if (difnew < .5*difold); epsred = .5*epsred; end;
    end;
    gv    = ztime (g, ipivot);
    gnorm = norm(gv, 'inf');
    ftest = 1 + abs(f);
    xnorm = norm(x,'inf');
    fprintf(1,'%4i   %4i   %4i   % .8e   %.1e\n', ...
        nit, nf, ncg, f, gnorm)
    %---------------------------------------------------------
    % test for convergence
    %---------------------------------------------------------
    [conv, flast, ipivot] = cnvtstm (alpha, pnorm, xnorm, ...
        difnew, ftest, gnorm, gtp, f, flast, g, ...
        ipivot, accrcy);
    conv_opt=conv;
    lamind=find(ipivot ~= 0 & ipivot ~= 2);
    toll = 10 * eps * (abs(low) + ones(size(x)));
    subl=find(x-low<toll);
    lambdaind1=intersect(lamind,subl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tolu = 10 * eps * (abs(up) + ones(size(x)));
    subu=find(up-x<tolu);
    lambdaindu=intersect(lamind,subu);
    lambdaind=struct('low',lambdaind1,'up',lambdaindu);
    
    if (nit>=maxiter); conv = 1; end;
    if (conv);
        if(conv_opt)
            disp('LMQNBC: termination 5');
        end
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    g = ztime (g, ipivot);
    %---------------------------------------------------------
    % modify data for LMQN preconditioner
    %---------------------------------------------------------
    if (~newcon);
        yk = g - oldg;
        sk = alpha*p;
        yksk = yk'*sk;
        ireset = (icycle == n-1 | difnew < epsred*(fkeep-f));
        if (~ireset);
            yrsr = yr'*sr;
            ireset = (yrsr <= 0);
        end;
        upd1 = (yksk <= 0);
    end;
    %---------------------------------------------------------
    % compute the search direction
    %---------------------------------------------------------
    argvec = [accrcy gnorm xnorm];
    [p, gtp, ncg1, d, eig_val{nit+1}] = ...
        modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
    ncg = ncg + ncg1;
    %---------------------------------------------------------
    % update LMQN preconditioner
    %---------------------------------------------------------
    if (~newcon);
        if (ireset);
            sr     = sk;
            yr     = yk;
            fkeep  = f;
            icycle = 1;
        else
            sr     = sr + sk;
            yr     = yr + yk;
            icycle = icycle + 1;
        end
    end;
end;
