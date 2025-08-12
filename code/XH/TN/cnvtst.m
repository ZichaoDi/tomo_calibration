function [conv, flast1, ipivot1] = cnvtst (alpha, pnorm, xnorm, ...
		 dif, ftest, gnorm, gtp, f, flast, g, ipivot, accrcy)
global lambdaind 
%---------------------------------------------------------
% test for convergence
%---------------------------------------------------------
% set up
%---------------------------------------------------------

toleps = sqrt(accrcy) + sqrt(eps);
rtleps = accrcy + eps;
conv   = 0;
imax   = 0;
ltest  = (flast - f <= -0.5*gtp);
%---------------------------------------------------------
% if anti-zigzag test satisfied, test multipliers;
% if appropriate, modify active set, only reduce active set
%---------------------------------------------------------
if (~ltest)
   ind   = find(ipivot ~= 0 & ipivot ~= 2);
   if (~isempty(ind));
      t = -ipivot(ind).*g(ind);
      [cmax, imax] = min(t);
      if (cmax >= 0); imax = 0; end;
   end;
end;
if (imax ~= 0)
   ipivot(ind(imax)) = 0;
   flast = f;
else
   stop1 = pnorm < toleps*(1 + xnorm);%alpha*pnorm < toleps*(1 + xnorm);
   stop2 = abs(dif) < rtleps*ftest;
   stop3 = gnorm < accrcy^(1/3)*ftest;
   stop4 = gnorm < .01*sqrt(accrcy)*ftest;
   stop5 = gnorm <1e-5;
   conv = ( stop1 ...
           & stop2        ...
           & stop3) ...
          | stop4;
          conv=stop5;
   if(conv)
       if(stop1 & stop2 & stop3)
           disp('tighter stoppting')
       elseif(stop5)
           disp('absolute gnorm')
       else
           disp('gnorm < eps * |f|')
       end
   end
end;
flast1  = flast;
ipivot1 = ipivot;
lambdaind   = find((ipivot ~= 2) & (ipivot.*g>0));
