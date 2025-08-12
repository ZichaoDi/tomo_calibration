function [x, f, g, nf1, alpha1, ierror, varargout] = ...
    lin_mg (p, x, f, alpha, g, sfun, low, up)
%---------------------------------------------------------
% line search: P(x+alpha*p)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ierror = 3;
maxit  = 5;
%%%%%%%%%%%%%%%%%%%#######################################################
if (alpha == 0); ierror = 0; maxit = 1; end;
alpha1 = alpha;
trialLength=4;
if(alpha1<=1)
    trialAlpha=linspace(1,alpha1,trialLength);
end
iproj=0;
q0=p'*g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
if(alpha1<=1)
    for trial=1:length(trialAlpha)
        xt = x + trialAlpha(trial).*p;
        [~,~, xt] = crash (xt, low, up);
        [ft, gt] = feval (sfun, xt);
        Armijo =ft<f+1e-4*trialAlpha(trial)*q0;
        if (Armijo );%& Wolfe);
            fprintf('Armijo and Wolfe satisfied, trial= %d\n',trial);
            ierror = 0;
            iproj  = 1;
            x   = xt;
            f   = ft;
            g   = gt;
            alpha1=trialAlpha(trial);
            itcnt=trial;
            break;
        end;
    end
    
end
if (alpha1 == 0); ierror = 0; maxit = 1; end;
if(iproj==0)
    for itcnt = 1:maxit;
        xt = x + alpha1.*p;
        [ft, gt] = feval (sfun, xt);
        if (ft < f);
            ierror = 0;
            x   = xt;
            f   = ft;
            g   = gt;
            break;
        end;
        alpha1 = alpha1 ./ 2;
    end;
    itcnt=itcnt+trialLength;
end
if (ierror == 3); alpha1 = 0; end;
nf1 = itcnt;

if (nargout == 7)
  dfdp = dot(gt, p);
  varargout{1} = dfdp;
end
