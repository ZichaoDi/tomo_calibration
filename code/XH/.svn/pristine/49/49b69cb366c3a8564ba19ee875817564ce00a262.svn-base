function [ipivot1, flast1] = modz (x, p, ipivot, low, up, flast, f, alpha)
%---------------------------------------------------------------------
% update the constraint matrix if a new constraint is encountered
%---------------------------------------------------------------------
indl = find(ipivot == 0 & p < 0);
if (length(indl) > 0);
   toll = 10 * eps * (abs(low(indl)) + ones(size(indl)));
   hitl = find(x(indl)-low(indl) <= toll);
   if (length(hitl)>0);
      flast = f;
      ipivot(indl(hitl)) = -1;
   end;
end;
%---------------------------------------------------------------------
indu = find(ipivot == 0 & p > 0);
if (length(indu) > 0);
   tolu = 10 * eps * (abs( up(indu)) + ones(size(indu)));
   hitu = find(up(indu)-x(indu)  <= tolu);
   if (length(hitu)>0 & norm(up)~=inf);
      flast = f;
      ipivot(indu(hitu)) = 1;
   end;
end;
%---------------------------------------------------------------------
flast1  = flast;
ipivot1 = ipivot;

