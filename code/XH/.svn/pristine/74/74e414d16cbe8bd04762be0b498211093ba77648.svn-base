function [xnew, fnew, gnew, nf1, alpha, ierror, varargout] = ...
	lin2 (p, x, f, alphamax, g, sfun)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
alpha  = alphamax;
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
%=======================================================
%===========test merit function=========================

% alp=linspace(0,2,20);
% for i=1:length(alp)
%     vt=x+alp(i)*p;
%     [F,G]=feval(sfun,vt);
%     [fm,gm]=merit_tn(vt,F,G);
%     meritf(i)=fm;
%     objf(i)=F;
% end
% figure(22);%subplot(5,2,it);
% plot(alp,meritf,'r*-',alp,objf,'b.-'); title('red: merit function; blue: original function');%pause;
% %==========================
for itcnt = 1:20;
   xt = x + alpha*p;
   [ft, gt] = feval (sfun, xt);
   if (ft < f);
      ierror = 0;
      xnew   = xt;
      fnew   = ft;
      gnew   = gt;
      fprintf('MG/Opt line search: alpha = %16.1e\n',alpha)
      break;
   end;
   alpha = alpha / 2;
end;

nf1 = itcnt;

if (nargout == 7)
  dfdp = dot(gt, p);
  varargout{1} = dfdp;
end
