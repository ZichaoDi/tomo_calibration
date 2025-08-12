function w = OPblurinpaint_pattern(u,A,S,tflag)

[m1,n1] = size(A);
[m2,n2] = size(S);
N = length(u); n = sqrt(N); 
if strcmpi(tflag,'size')
    w(1) = m1*m2;
    w(2) = n1*n2;
elseif strcmpi(tflag,'notransp')
    w = A*u;
    w = reshape(w, n, n);
    w = S.*w;
    w = w(:);
else
    u = reshape(u, n, n);
    w = S.*u;
    w = A'*w(:);
end