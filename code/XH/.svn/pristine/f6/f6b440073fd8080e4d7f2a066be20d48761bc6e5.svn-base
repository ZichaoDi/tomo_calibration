function [xnew, fnew, gnew, nf1, ierror, alpha1, varargout] = ...
    lin1 (p, x, f, alpha, g, sfun)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
maxit  = 15;
%%%%%%%%%%%%%%%%%%%#####fix step length as 1 ##############################
%   alpha=1;
%%%%%%%%%%%%%%%%%%%#######################################################
% disp('### In LINE SEARCH: lin1 ###')
% fprintf(' g''p = %e\n',g'*p);
if (alpha == 0); ierror = 0; maxit = 1; end;
alpha1 = alpha;
q0=p'*g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
for itcnt = 1:maxit;
    xt = x + alpha1.*p;
    %%%%%%%%%%%#############################################################
    [ft, gt] = feval (sfun, xt);
    % Armijo =ft<f+1e-4*alpha1*q0;
    % Wolfe = abs(p'*gt)<0.25*abs(q0);
    if (ft < f);
%     if(Armijo & Wolfe)
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        gnew   = gt;
        break;
    end;
    alpha1 = alpha1 ./ 2;
end;
if (ierror == 3); alpha1 = 0; end;
nf1 = itcnt;

if (nargout == 7)
    dfdp = dot(gt, p);
    varargout{1} = dfdp;
end
