function td = initpc (d, upd1, ireset)
%---------------------------------------------------------
% initialize the diagonal preconditioner
%---------------------------------------------------------
global hyk ykhyk hyr yksr ykhyr yrhyr sk yk sr yr ...
       hg gsk yksk yrsr
%---------------------------------------------------------
if (upd1)
   td = d;
else
   if (ireset)
      hyk  = d.*sk;
      sds  = sk'*hyk;
      if (hyk == 0);
        disp('INITPC: hyk = 0')
        nd = norm(d)
        nsk = norm(sk)
      end
      if (sds == 0);
        disp('INITPC: sds = 0')
        nsk = norm(sk)
        nhyk = norm(hyk)
      end;
      td   = d - d.*d.*sk.*sk/sds + yk.*yk/yksk;
   else
      hyk  = d.*sr;
      sds  = sr'*hyk;
      srds = sk'*hyk;
      yrsk = yr'*sk;
      hyk  = d.*sk - hyk*srds/sds + yr*yrsk/yrsr;
      td   = d - d.*d.*sr.*sr/sds+yr.*yr/yrsr;
      sds  = sk'*hyk;
      td   = td - hyk.*hyk/sds + yk.*yk/yksk;
   end
end;
