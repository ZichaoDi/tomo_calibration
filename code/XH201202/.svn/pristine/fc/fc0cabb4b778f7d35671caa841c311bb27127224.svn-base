function [spe]= stpmax (stepmx, pe, x, p, ipivot, low, up)
%------------------------------------------------
% compute the maximum allowable step length
% (spe is the standard (unconstrained) max step)
%------------------------------------------------
spe  = stepmx / pe;
al   = spe;
au   = spe;
%------------------------------------------------
indl = find(ipivot==0 & p < 0);% 
if (~isempty(indl));
   tl   = low(indl) - x(indl);
   [al,indsub]  = min(tl./p(indl));

end;
%------------------------------------------------
indu = find( ipivot==0 & p > 0);%
if (~isempty(indu) );
   tu   = up (indu) - x(indu);
   au   = min(tu./p(indu));
end;
%------------------------------------------------
spe  = min([spe al au]);

