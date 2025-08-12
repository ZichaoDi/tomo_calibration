function [xnew, fnew, gnew, nf1, ierror, alpha] = ...
    lin_Cauchy (p, x, f, g, sfun, alphaDis)
%------------------------------------------------------------------------
% line search along the piecewise linear path starting from Cauchy point
%------------------------------------------------------------------------
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
maxit  = 15;
nf1=0;
%%%%%%%%%%%%%%%%%%%#######################################################
if (alphaDis == 0); ierror = 0; maxit = 1; end;
jmax=length(alphaDis);
%---------------------------------------------------------
% line search
%---------------------------------------------------------
for itcnt = jmax:-1:2;
    d=p(:,itcnt-1);
    q0=d'*g;
    %%%%%%%%%%%#############################################################
    if(itcnt==1)
        alpha_min=0;
    else
        alpha_min=alphaDis(itcnt-1);
    end
    for it=1:maxit
        xt = x + alphaDis(itcnt).*d;
        [ft, gt] = feval (sfun, xt);
        Armijo =ft<f+1e-4*alphaDis(itcnt)*q0;
        Wolfe = abs(d'*gt)<0.25*abs(q0);
        if(Armijo)
            ierror = 0;
            xnew   = xt;
            fnew   = ft;
            gnew   = gt;
            alpha=alphaDis(itcnt);
            break;
            
        elseif(alphaDis(itcnt)<=alpha_min)
            break;
        end;
        alphaDis(itcnt) = alphaDis(itcnt)/2;
    end
        nf1 = nf1+it;
end;
if (ierror == 3); alpha = 0; end;
