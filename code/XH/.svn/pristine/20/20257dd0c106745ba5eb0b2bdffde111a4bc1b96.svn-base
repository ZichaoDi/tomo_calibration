function ipivot=pivot_naive(x,g,low,up)
tol=eps;
ipivot = zeros(size(x));
indl=find(abs(x-low)<=tol & g>0);
if(~isempty(indl))
    ipivot(indl) =  -1;
end
indu=find(abs(x-up)<=tol & g<0);
if(~isempty(indu))
    ipivot(indu) =  1;
end
