function [xstar, f, g, ierror] = ...
    lmqnbc (x, sfun, low, up, maxit, maxfun, stepmx, accrcy)
%---------------------------------------------------------
% This is a bounds-constrained truncated-newton method.
% The truncated-newton method is preconditioned by a
% limited-memory quasi-newton method (computed by
% this routine) with a diagonal scaling (routine ndia3).
% For further details, see routine tnbc.
%---------------------------------------------------------
global shift
global sk yk sr yr yksk yrsr
global NF N current_n x_iter fiter itertest ErrIter
global ptest gv ipivot nit
global n_delta N_delta i_cauchy W0  m NumElement
global maxiter err0 Joint 
global f_xrf f_xtm
%---------------------------------------------------------
% check that initial x is feasible and that the bounds
% are consistent
%---------------------------------------------------------
[f, g] = feval (sfun, x);
oldf   = f;
fiter=[];
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
if(Joint==1)
fprintf(1,'  it     nf     cg       f         f_xrf        f_xtm       |g|       alpha      error\n');
else
fprintf(1,'  it     nf     cg       f          |g|        alpha          error\n');   
end

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
if(Joint==1)
    [f,g, f_xrf, f_xtm] = feval (sfun, x);
else
    [f,g] = feval (sfun, x);
end
g0=g;
nf     = 1;
nit    = 0;
nitOld=nit;
flast  = f;
itertest=[];
ErrIter=[];
x_iter=[];
x_iter(:,1)=x;
itertest(1)=nf+ncg;
ErrIter(1)=err0;
ind = find((ipivot ~= 2) & (ipivot.*g>0));
if (~isempty(ind));
    ipivot(ind) = zeros(length(ind),1);
end;
ipivotOld=ipivot;
g = ztime (g, ipivot);
gnorm = norm(g,'inf');
if(Joint==1)
fprintf(1,'%4i   %4i   %4i   %.3e   %.3e    %.3e    %.1e    %.1e   %.3e\n', ...
    nit, nf, ncg, f,  f_xrf, f_xtm,gnorm, 1, err0);
else
    fprintf(1,'%4i   %4i   %4i   % .3e   %.1e     %.1e      %.3e\n', ...
        nit, nf, ncg, f, gnorm, 1, err0);
end
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
i_cauchy=0;
d      = ones(n,1);
%---------------------------------------------------------
% ..........main iterative loop..........
%---------------------------------------------------------
% compute the search direction
%---------------------------------------------------------
argvec = [accrcy gnorm xnorm];
%%%%%%%%%%%%%%%%%%%%%##############################
%%%%%%%%%%%%%%%%%%%%%##############################
[p, gtp, ncg1, d, eig_val] = ...
    modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
ncg = ncg + ncg1;
% figure,
while (~conv);
    oldg = g;
    pnorm = norm(p, 'inf');
    oldf = f;
    fiter(nit+2)=oldf;
    itertest(nit+2)=nf+ncg;
    % ang=acos(dot(g,p)/(norm(g)*norm(p)));
    % ang*180/pi
    % ang=acos(dot(g(1:n_delta),p(1:n_delta))/(norm(g(1:n_delta))*norm(p(1:n_delta))));
    % ang*180/pi
    
    %---------------------------------------------------------
    % line search
    %---------------------------------------------------------
    pe = pnorm + eps;
    spe = stpmax (stepmx, pe, x, p, ipivot, low, up);
    %---------------------------------------------------------%    
    % update active set, if appropriate
    %---------------------------------------------------------
    alpha = step1 (f, gtp, spe);
    alpha0 = alpha;
    PieceLinear=1;
    newcon = 0;
    if(PieceLinear)
        if(Joint==1)
        [x_new, f_new, g_new, nf1, ierror, alpha,ipivot,newcon,flast,f_xrf,f_xtm] = lin_proj (p, x, f, g, alpha0, sfun, low, up,ipivot,newcon,flast);
        else
        [x_new, f_new, g_new, nf1, ierror, alpha,ipivot,newcon,flast] = lin_proj (p, x, f, g, alpha0, sfun, low, up,ipivot,newcon,flast); 
        end
    else
        [x_new, f_new, g_new, nf1, ierror, alpha] = lin1 (p, x, f, alpha0, g, sfun);
    end

    if (alpha <= 0 & alpha0 ~= 0 | ierror == 3);
        fprintf('Error in Line Search\n');
        fprintf('    ierror = %3i\n',    ierror);
        fprintf('    alpha  = %12.6f\n', alpha);
        fprintf('    alpha0 = %12.6f\n', alpha0);
        fprintf('    g''p    = %12.4e\n', gtp);
        %############################
        fprintf('    |g|     = %12.4e\n', norm(g));
        fprintf('    |p|     = %12.4e\n', norm(p));
    end;
    %#######################
    x   = x_new;
    f   = f_new;
    g   = g_new;
    g0  = g;
    %#######################
    nf  = nf  + nf1;
    nit = nit +   1;
%---------------------------------------------------------
% update active set, if appropriate
%---------------------------------------------------------
    newcon = 0;
    if (abs(alpha-spe) <= 10*eps);% | alpha==1
        % disp('update ipivot due to tiny step length')
        newcon = 1;
        ierror = 0;
        [ipivot, flast] = modz (x, p, ipivot, low, up, flast, f);
    end;
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
    x_iter(:,nit+1)=x;
    %--------------------------------- Error
    ErrIter(nit+1)=norm(x-W0);
    if(Joint==1)
        fprintf(1,'%4i   %4i   %4i   %.3e   %.3e    %.3e    %.1e    %.1e   %.3e\n', ...
        nit, nf, ncg, f, f_xrf, f_xtm, gnorm, alpha, norm(x-W0));
    else
        fprintf(1,'%4i   %4i   %4i   % .6e   %.2e     %.1e      %.3e\n', ...
        nit, nf, ncg, f, gnorm, alpha, norm(x-W0));
    end
    %---------------------------------------------------------
    % test for convergence
    %---------------------------------------------------------
    [conv, flast, ipivot] = cnvtst (alpha, pnorm, xnorm, ...
        difnew, ftest, gnorm, gtp, f, flast, g, ...
        ipivot, accrcy);
    if(nit>=maxiter)
        conv = 1;
    end;

    plotAS=0;
    if((mod(nit,1)==0 | conv | nit==1) & plotAS==1) %
        figure(10);
        subplot(2,2,1);imagesc(reshape(x,m(1),m(2)));
        colorbar;title('x');
        subplot(2,2,2);imagesc(reshape(ipivot,m(1),m(2)));
        colorbar;title('ipivot');
        subplot(2,2,3);imagesc(reshape(sign(g),m(1),m(2)));
        colorbar;title('gradient');
        subplot(2,2,4);imagesc(reshape(sign(p),m(1),m(2)));
        colorbar;title([num2str(nit) ,'direction']);drawnow;

        % addAS=length(find(-ipivot+ipivotOld==1));
        % dropAS=length(find(-ipivot+ipivotOld==-1));
        % figure(91);subplot(2,1,1),
        % plot(1:length(x),ipivot,'r.-',1:length(x),x,'bo-',1:length(x),g,'gs');
        % ipivotOld=ipivot;
        % legend('active set','current variable','reduced gradient');
        % subplot(2,1,2);qpPlotAset(ipivot,nit,length(x),[addAS,-dropAS],nitOld);
        % nitOld=nit;
        % if(conv)
        %     figure(90);clims=[-1 1];
        %     for inum=1:NumElement
        %         ASpar=4;
        %         subplot(NumElement,ASpar,ASpar*(inum-1)+1)
        %         AS=find(x-0<=10 * eps);
        %         xmap=zeros(size(x));xmap(AS)=-1;xmap=reshape(xmap(n_delta+1:end),m(1),m(2),NumElement);
        %         imagesc(xmap(:,:,inum),clims);ylabel(['Element',num2str(inum)],'fontsize',12);
        %         if(inum==1);title('Indicator of x','fontsize',12);end
        %         gmap=g;
        %         gmap(AS)=-1*sign(g(AS));gmap=reshape(gmap(n_delta+1:end),m(1),m(2),NumElement);
        %         subplot(NumElement,ASpar,ASpar*(inum-1)+2)
        %         imagesc(gmap(:,:,inum),clims);
        %         if(inum==1);title('Indicator of gradient','fontsize',12);end
        %             ipmap=reshape(ipivot(n_delta+1:end),m(1),m(2),NumElement);
        %         subplot(NumElement,ASpar,ASpar*(inum-1)+3)
        %         imagesc(ipmap(:,:,ikum),clims);
        %         if(inum==1);title('Indicator of active set','fontsize',12);end
        %         
        %         pmap=p;
        %         pmap(AS)=sign(p(AS));pmap=reshape(pmap(n_delta+1:end),m(1),m(2),NumElement);
        %         subplot(NumElement,ASpar,ASpar*(inum-1)+4)
        %         imagesc(pmap(:,:,inum),clims);
        %         if(inum==1);title('search direction','fontsize',12);end
        %     end
        %     hp4 = get(subplot(NumElement,ASpar,ASpar*NumElement),'Position');
        %     colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1])
        % end
        drawnow;
        pause;
    end
    %------------------------------------------------------
    if (conv);
        disp('LMQNBC: termination 5')
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    
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
    
    [p, gtp, ncg1, d,~,cgstop] = ...
        modlnp (d, x, g, 10, upd1, ireset, bounds, ipivot, argvec, sfun);
    ptest=p;
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
