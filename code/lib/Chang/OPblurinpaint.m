function w = OPblurinpaint(u,A,S,tflag)

[m1,n1] = size(A);
[m2,n2] = size(S);
if strcmpi(tflag,'size')
    w(1) = m1*m2;
    w(2) = n1*n2;
elseif strcmpi(tflag,'notransp')
    w = A*u;
    w = S*w;
else
    w = S'*u;
    w = A'*w;
end