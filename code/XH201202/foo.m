function gh=foo(fctn,x0);
global numThetan testind
z0 = x0;
z = z0(:);
n = length(z);
h =1;
testind=0;
[f,g] =feval(fctn,z0);
diff_ind = 0; % Use 1 for central differencing, 0 for forward differencing
gh=zeros(size(g));
if (diff_ind == 0);
%     disp('Forward Differencing')
%     fprintf(' i                g                gh         error\n')
    for i= 1:500:n;
        testind=i;
        zh = z;
        zh(i) = zh(i) + h;
        fh = feval(fctn,zh);
        gh(i,:) = (fh(:)-f(:))./h;
        eee = norm(g(i,:)-gh(i,:))/(1+norm(g(i,:)));
        fprintf('%3i     %8.8e        %8.8e     %8.8e\n',i,norm(g(i,:),inf),norm(gh(i,:),inf),eee);
    end;
     err = norm(g-gh)/(1+norm(g));
else
    disp('Central Differencing')
    fprintf(' i               fh1                  fh2                  g                       gh              error\n')
    for i= 1:n;
        testind=i;
        zh = z;
        zh(i) = zh(i) + h;
        fh1 = feval(fctn,zh);
        zh = z;
        zh(i) = zh(i) - h;
        fh2 = feval(fctn,zh);
        gh(i,:) = (fh1-fh2)./(2*h);
        eee = norm(g(i,:)-gh(i,:))/(1+norm(g(i,:)));
        fprintf('%3i     %18.8e  %18.8e  %18.8e  %18.8e     %8.1e\n',i,norm(fh1),norm(fh2),norm(g(i,:)),norm(gh(i,:)),eee);
        %     ee2=norm(fh1-2*f+fh2,'inf')/h^2
    end;
    err = norm(g-gh,'inf')/(1+norm(g,'inf'))
end;
