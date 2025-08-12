function y = msolve (g, upd1, ireset, first, d)
%---------------------------------------------------------
% This routine acts as a preconditioning step for the
% linear conjugate-gradient routine.  It is also the
% method of computing the search direction from the
% gradient for the non-linear conjugate-gradient code.
% It represents a two-step self-scaled bfgs formula.
%---------------------------------------------------------
global hyk ykhyk hyr yksr ykhyr yrhyr sk yk sr yr ...
       hg gsk yksk yrsr
%---------------------------------------------------------
if (upd1)
   y = g./d;
else
   gsk = g'*sk;
   if (ireset)
      hg = g./d;
      if (first);
	 hyk   = yk./d;
         ykhyk = yk'*hyk;
      end;
      ghyk = g'*hyk;
      y    = ssbfgs(sk,hg,hyk,yksk,ykhyk,gsk,ghyk);
   else
      hg = g./d;
      if (first);
         hyk   = yk./d;
         hyr   = yr./d;
         yksr  = yk'*sr;
         ykhyr = yk'*hyr;
      end;
      gsr = g'*sr;
      ghyr = g'*hyr;
      if (first);
         yrhyr = yr'*hyr;
      end;
      hg = ssbfgs(sr,hg,hyr,yrsr,yrhyr,gsr,ghyr);
      if (first);
	  hyk = ssbfgs(sr,hyk,hyr,yrsr,yrhyr,yksr,ykhyr);
      end;
      ykhyk = hyk'*yk;
      ghyk  = hyk'*g;
      y     = ssbfgs(sk,hg,hyk,yksk,ykhyk,gsk,ghyk);
   end
end;
