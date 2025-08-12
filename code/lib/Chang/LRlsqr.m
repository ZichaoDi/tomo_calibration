function [X, relerr,relres,rkX, lambdahis,U,V,M,T,Z] = LRlsqr(A, b, x0,xtrue, parameters)

% input:
% iftrun - 1 for computing truncated SVD of V and X, 2 for doing singular 
%           value shrinkage (as in SVT)
% tau - shrinkage parameter (as in SVT)
% reg - regularization parameter for hybrid FLSQR; 'discrep' for discrepancy 
%           principle based parameter selection, or other numerical values
% eta, nl - used for discrepancy principle; eta is slightly larger than 1,
%           nl is noise level



k = parameters.maxIt;

kappaB = parameters.kappaB;
kappa = parameters.kappa;
iftrun = parameters.iftrun;
% thr = parameters.thr;

reg = parameters.reg;
eta = parameters.eta;
nl = parameters.nl;



tau = parameters.tau;

m = length(b);
n = length(Atransp_times_vec(A,ones(m,1)));
normx = norm(xtrue);
normb = norm(b);

r0 = b - A_times_vec(A,x0);
beta1 = norm(r0); 
u1 = r0/beta1;

Z = zeros(n,k);
M = zeros(k+1,k);
T = zeros(k+1);
U = zeros(m,k+1); U(:,1) = u1;
V = zeros(n,k+1);

X = zeros(n,k);
relerr = zeros(k,1);
relres = relerr;
rkX = relerr;

nim = sqrt(n);

etaeps = eta*nl;

lambdahis = zeros(k,1);

if strcmp(reg,'discrep')
    lambda = 1;
elseif isnumeric(reg)
    lambda = reg;
end




for i = 1:k
    
    w = Atransp_times_vec(A,U(:,i));
    
    for j = 1:i-1
        T(j,i) = w'*V(:,j);
        w = w - T(j,i)*V(:,j);
    end
    
    T(i,i) = norm(w);
    V(:,i) = w/T(i,i);
    
    Z(:,i) = V(:,i);
    
    if iftrun == 1
        z = truncfcn(V(:,i),sqrt(n),kappaB);
        Z(:,i) = z;
%         z = reshape(V(:,i),nim,nim);
%         [UZ,SZ,VZ]=svd(z);
%         sz = diag(SZ);
%         diffsZ = sz(1:end-1)-sz(2:end); 
%         diffsZr = diffsZ./diffsZ(1);
%         trunc = find(diffsZr<thr, 1)
%         
%         sz(trunc+1:end) = 0;
%         z = UZ*diag(sz)*VZ';
%         Z(:,i) = z(:);
    elseif iftrun == 2
        z = softthr(V(:,i),sqrt(n),tau);
        Z(:,i) = z;
    end
    
    w = A_times_vec(A,Z(:,i));
        
    for j = 1:i
        M(j,i) = w'*U(:,j);
        w = w - M(j,i)*U(:,j);
    end
    
    M(i+1,i) = norm(w);
    U(:,i+1) = w/M(i+1,i);
    
    d = zeros(i+1,1);
    d(1) = beta1;
    
%     y = [M(1:i+1,1:i);lambda*eye(i)]\[d;zeros(i,1)];
    
    [UM,SM,VM] = svd(M(1:i+1,1:i));
    lsqr_res = abs((UM(:,i+1)'*d))/normb;
    
    c = UM'*d;
    
    if i == 1
        SM = SM(1,1);
    else
        SM = diag(SM);
    end
    
    
    if iftrun == 0 && strcmp(reg,'opt')
        xhat = VM'*(V(:,1:i)'*xtrue);
        sB = SM;
        myfun = @(lambda) norm(((conj(sB)./(abs(sB).^2+lambda)).*c(1:end-1) - xhat));
        options.TolX = eps;
        lambda = fminbnd(myfun, 0, sB(1),options);
        lambdahis(i) = lambda; 
    end

        
    
    Filt = SM.^2 + lambda;
    y = VM*((SM.*c(1:i))./Filt);
    
    %%% NOTE: the residual is computed before performing the last
    %%% truncation (of the solution)!!!
    relres(i) = norm(M(1:i+1,1:i)*y - d)/normb;
    

    if strcmp(reg,'discrep')
        lambda = abs((etaeps-lsqr_res)/(relres(i)-lsqr_res))*lambda;
        lambdahis(i+1) = lambda; 
    end    
    

    x = x0 + Z(:,1:i)*y;
    X(:,i) = x;
    
    if iftrun == 1
        X(:,i) = truncfcn(X(:,i),sqrt(n),kappa);
    elseif iftrun == 2
        X(:,i) = softthr(X(:,i),sqrt(n),tau);
    end
    
    relerr(i) = norm(xtrue - X(:,i))/normx;   
    
    Xtemp = reshape(X(:,i), sqrt(n), sqrt(n));
    sX = svd(Xtemp);
    rkX(i) = length(sX(sX>1e-10));
    
end        
        