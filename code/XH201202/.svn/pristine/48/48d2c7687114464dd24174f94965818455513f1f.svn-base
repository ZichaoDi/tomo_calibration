function [xstar, f, g, ierror, eig_val] = ...
    lmqnbcm1 (x, sfun, low, up, maxiter, maxit, maxfun, stepmx, accrcy)
%---------------------------------------------------------
% This is a bounds-constrained truncated-newton method.
% The truncated-newton method is preconditioned by a
% limited-memory quasi-newton method (computed by
% this routine) with a diagonal scaling (routine ndia3).
% For further details, see routine tnbc.
%---------------------------------------------------------
global sk yk sr yr yksk yrsr
global NF N current_n  fiter itertest
global gv ipivot nit
global i_cauchy W0  m NumElement
global Joint 
%---------------------------------------------------------
% check that initial x is feasible and that the bounds
% are consistent
%---------------------------------------------------------
[f, g] = feval (sfun, x);
oldf   = f;
fiter(1)=oldf;
[ipivot, ierror, x] = crash(x, low, up);
if (ierror ~= 0);
    disp('LMQNBC: termination 0')
    fprintf(1,'LMQNBC: terminating (no feasible point)');
    f = 0;
    g = zeros(size(x));
    xstar = x;
    return;
end;
%---------------------------------------------------------
% initialize variables, parameters, and constants
%---------------------------------------------------------
fprintf(1,'  it     nf     cg           f            |g|      alpha        error\n');   

nind   = find(N==current_n);
upd1   = 1;
ncg    = 0;
conv   = 0;
xnorm  = norm(x,'inf');
ierror = 0;
if (stepmx < sqrt(accrcy) | maxfun < 1);
    disp('LMQNBC: termination 1')
    ierror = -1;
    xstar = x;
    return;
end;
%---------------------------------------------------------
% compute initial function value and related information
%---------------------------------------------------------
[f, g] = feval (sfun, x);
nf     = 1;
nit    = 0;
nitOld=nit;
flast  = f;
itertest(1)=nf+ncg;
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
ipivotOld=ipivot;
gnorm = norm(g,'inf');
fprintf(1,'%4i   %4i   %4i   % .8e   %.1e     %.1e      %.3e\n', ...
    nit, nf, ncg, f, gnorm, 1, norm(x-W0));
%---------------------------------------------------------
% check if the initial point is a local minimum.
%---------------------------------------------------------
ftest = 1 + abs(f);
if (gnorm < .01*sqrt(eps)*ftest);
    disp('LMQNBC: termination 2')
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
eig_val{1} = [];
%---------------------------------------------------------
% ..........main iterative loop..........
%---------------------------------------------------------
% compute the search direction
%---------------------------------------------------------
argvec = [accrcy gnorm xnorm];
%%%%%%%%%%%%%%%%%%%%%##############################
[p, gtp, ncg1, d, eig_val{nit+1}] = ...
    modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
ncg = ncg + ncg1;
while (~conv);
    oldg = g;
    pnorm = norm(p, 'inf');
    oldf = f;
    fiter(nit+2)=oldf;
    itertest(nit+2)=nf+ncg;
    
    %---------------------------------------------------------
    % line search
    %---------------------------------------------------------
    pe = pnorm + eps;
    [spe] = stpmax (stepmx, pe, x, p, ipivot, low, up);
    alpha = step1 (f, gtp, spe);
    alpha0 = alpha;
    PieceLinear=1;
    newcon = 0;
    if(PieceLinear)
        [x_new, f_new, g_new, nf1, ierror, alpha,ipivot,newcon,flast] = lin_proj (p, x, f, g, alpha0, sfun, low, up,ipivot,flast,newcon); 
    else
        [x_new, f_new, g_new, nf1, ierror, alpha] = lin1 (p, x, f, alpha0, g, sfun);
    end
    
    Cauchy=0;
    %---------------------------------------------------------%    
    if (alpha == 0 & alpha0 ~= 0 | ierror == 3);
        fprintf('Error in Line Search\n');
        fprintf('    ierror = %3i\n',    ierror);
        fprintf('    alpha  = %12.6f\n', alpha);
        fprintf('    alpha0 = %12.6f\n', alpha0);
        fprintf('    g''p    = %12.4e\n', gtp);
        %############################
        fprintf('    |g|     = %12.4e\n', norm(g));
        fprintf('    |p|     = %12.4e\n', norm(p));
        %         disp('Hit any key to continue')
        %         pause;
    end;
    %#######################
    x   = x_new;
    f   = f_new;
    g   = g_new;
    g0=g;
    %#######################
    nf  = nf  + nf1;
    nit = nit +   1;
    ASchange(nit)=norm(-ipivot+ipivotOld,1);
    if (ierror == 3);
        disp('LMQNBC: termination 3')
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
        disp('LMQNBC: termination 4')
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
    %--------------------------------- Error
    ErrIter(nit)=norm(x-W0);
    fprintf(1,'%4i   %4i   %4i   % .8e   %.1e     %.1e      %.3e\n', ...
        nit, nf, ncg, f, gnorm, alpha, norm(x-W0));
    %---------------------------------------------------------
    % test for convergence
    %---------------------------------------------------------
    [conv, flast, ipivot] = cnvtst (alpha, pnorm, xnorm, ...
        difnew, ftest, gnorm, gtp, f, flast, g, ...
        ipivot, accrcy);
    if(nit>=maxiter)
        conv = 1;
    end;
    %------------------------------- Active set plot
    plotAS=0;
    if((mod(nit,10)==0 | conv | nit==1) & Joint~=-1 & plotAS==1) %
        %     if(conv& Joint~=-1 & plotAS)
        load BindInd
        figure(10);plot(1:length(x),g,'r.-',1:length(BindInd),g(BindInd),'ro',1:length(x),ipivot,'g.-',1:length(x),p,'m*');
        addAS=length(find(-ipivot+ipivotOld==1));
        dropAS=length(find(-ipivot+ipivotOld==-1));
        figure(19);plot(1:nit,ASchange,'r.-');
        figure(91);subplot(2,1,1),
        plot(1:length(x),ipivot,'r.-',1:length(x),x,'bo-',1:length(x),g,'gs');
        ipivotOld=ipivot;
        legend('active set','current variable','reduced gradient');
        subplot(2,1,2);qpPlotAset(ipivot,nit,length(x),[addAS,-dropAS],nitOld);
        nitOld=nit;
        if(conv)
            figure(90);clims=[-1 1];
            for inum=1:NumElement
                ASpar=4;
                subplot(NumElement,ASpar,ASpar*(inum-1)+1)
                AS=find(x-0<=10 * eps);
                xmap=zeros(size(x));xmap(AS)=-1;xmap=reshape(xmap,m(1),m(2),NumElement);
                imagesc(xmap(:,:,inum),clims);ylabel(['Element',num2str(inum)],'fontsize',12);
                if(inum==1);title('Indicator of x','fontsize',12);end
                gmap=g;
                gmap(AS)=-1*sign(g(AS));gmap=reshape(gmap,m(1),m(2),NumElement);
                subplot(NumElement,ASpar,ASpar*(inum-1)+2)
                imagesc(gmap(:,:,inum),clims);
                if(inum==1);title('Indicator of gradient','fontsize',12);end
                ipmap=reshape(ipivot,m(1),m(2),NumElement);
                subplot(NumElement,ASpar,ASpar*(inum-1)+3)
                imagesc(ipmap(:,:,inum),clims);
                if(inum==1);title('Indicator of active set','fontsize',12);end
                
                pmap=p;
                pmap(AS)=sign(p(AS));pmap=reshape(pmap,m(1),m(2),NumElement);
                subplot(NumElement,ASpar,ASpar*(inum-1)+4)
                imagesc(pmap(:,:,inum),clims);
                if(inum==1);title('search direction','fontsize',12);end
            end
            hp4 = get(subplot(NumElement,ASpar,ASpar*NumElement),'Position');
            colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1])
        end
        drawnow;
    end
    %------------------------------------------------------
    if (conv);
%         disp('LMQNBC: termination 5')
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    
    %%%######################## Update active set using Cauchy point
    if(Cauchy)
        [~,i_cauchy,alpha_cp,p_cp]=Cauchy_Point(x,g0, accrcy, xnorm, sfun,low,up);
        if(i_cauchy~=0)
            [x, f, g, nf1, ierror] = lin_Cauchy (p_cp, x, f, g0, sfun, alpha_cp);
            if(ierror==3)
                disp('descent not found from Cauchy point')
            else
                nf = nf+nf1;
                [ipivot] = crash2 (x,g0, low, up);
            end
        else
            disp('Cauchy point not found')
        end
    end
    %%%===========================================================
    
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
        modlnp (d, x, g, 10, upd1, ireset, bounds, ipivot, argvec, sfun);
    ncg = ncg + ncg1;
    % %---------------------------------------------------------
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
