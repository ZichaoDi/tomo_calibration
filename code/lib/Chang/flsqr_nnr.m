function [X, relerr,relres,objhis,lambdahis,alphahis,U,V,M,T,Z] = flsqr_nnr(A, b, x0,xtrue, parameters)

p = parameters.p;

k = parameters.maxIt;
svdbasis = parameters.svdbasis;
% thr = parameters.thr;

reg = parameters.reg;
eta = parameters.eta;
nl = parameters.nl;
gamma = parameters.gamma;

m = length(b);
n = length(Atransp_times_vec(A,ones(m,1)));
normxtrue = norm(xtrue);
normb = norm(b);

nim = sqrt(n);

% if m==n
%     B = reshape(b,nim,nim);
%     [UB,SB,VB] = svd(B);
%     UX = UB; VX = VB; SX = SB;% in the case x0 = 0;
% else
%     UX = eye(nim); VX = eye(nim); 
% end




r0 = b - A_times_vec(A, x0);
beta1 = norm(r0); 
u1 = r0/beta1;

UX = eye(nim); VX = eye(nim);
SX = eye(nim); 
sx = diag(SX);

% if strcmp(svdbasis,'t') % truncate small singular values
%     diffsX = sx(1:end-1)-sx(2:end); 
%     diffsXr = diffsX./diffsX(1);
%     trunc = find(diffsXr<thr, 1);
%     sx(trunc+1:end) = 0;
%     precX = kron(ones(n,1),sx);
%     precX = precX.^((2-p));
% else
if strcmp(svdbasis,'x')
    precX = kron(ones(n,1),sx.^2);
    precX = precX + gamma;
    precX = precX.^((2-p)/2);
end

Z = zeros(n,k);
M = zeros(k+1,k);
T = zeros(k+1);
U = zeros(m,k+1); U(:,1) = u1;
V = zeros(n,k+1);

X = zeros(n,k);
relerr = zeros(k,1);
relres = relerr;
objhis = relerr;

lambdahis = zeros(k+1,1);
alphahis = zeros(k,1);

if strcmp(reg,'discrep')
    lambda = parameters.RegParam0;
elseif isnumeric(reg)
    lambda = reg;
end

etaeps = eta*nl;


for i = 1:k
    
    w = Atransp_times_vec(A,U(:,i));
          
    for j = 1:i-1
        T(j,i) = w'*V(:,j);
        w = w - T(j,i)*V(:,j);
    end
    
    T(i,i) = norm(w);
    V(:,i) = w/T(i,i);

    
%%%%
      vtemp = V(:,i);
      vtemp = reshape(vtemp, nim, nim);
      
      if strcmp(svdbasis,'v')
        [UV, SV, VV] = svd(vtemp);
%         sv = diag(SV);
% %         diffsV = sv(1:end-1)-sv(2:end); 
% %         diffsVr = diffsV./diffsV(1);
% %         trunc = find(diffsVr<thr, 1);
% %         sv(trunc+1:end) = 0;

%         precX = kron(ones(nim,1), sv);
%         precX = precX.^((2-p));
%         vtemp = UV'*(vtemp*VV);
%         vtemp = precX.*vtemp(:);
        precX = SV.^(3-p);
        vtemp = precX(:);
        vtemp = reshape(vtemp, nim, nim); 
        vtemp = UV*(vtemp*VV');
      else
        if i > 1
        vtemp = UX'*(vtemp*VX);
        vtemp = precX.*vtemp(:);
        vtemp = reshape(vtemp, nim, nim); 
        vtemp = UX*(vtemp*VX');
        end
      end
      
      
      Z(:,i) = vtemp(:);
      [Qtemp, Rtemp] = qr(Z(:,1:i), 0);
%%%%
    
    w = A_times_vec(A,Z(:,i));
    

    for j = 1:i
        M(j,i) = w'*U(:,j);
        w = w - M(j,i)*U(:,j);
    end
    
    M(i+1,i) = norm(w);
    U(:,i+1) = w/M(i+1,i);

    
    d = zeros(i+1,1);
    d(1) = beta1;
    
    
    [UM,SM,VM] = svd(M(1:i+1,1:i));
    lsqr_res = abs((UM(:,i+1)'*d))/normb;
    
%     c = UM'*d;
%     
%     if i == 1
%         SM = SM(1,1);
%     else
%         SM = diag(SM);
%     end
%     
%     Filt = SM.^2 + lambda; Filt0 = SM.^2;
%     y = VM*((SM.*c(1:i))./Filt);
%     
%     y0 = VM*((SM.*c(1:i))./Filt0);
%     alphahis(i) = norm(M(1:i+1,1:i)*y0 - d)/normb;

    y = [M(1:i+1,1:i); sqrt(lambda)*Rtemp]\[d; zeros(i,1)];
    y0 = M(1:i+1,1:i)\d;
    alphahis(i) = norm(M(1:i+1,1:i)*y0 - d)/normb;
    
    resi = norm(M(1:i+1,1:i)*y - d);
    relres(i) = resi/normb;
    
    
    x = x0 + Z(:,1:i)*y;
    X(:,i) = x;
    
    XX = reshape(x,nim,nim);
    sx = svd(XX);
    objhis(i) = resi^2 + lambda*sum(sx);
    
    
    
    
    if strcmp(reg,'discrep')
%         etaeps
%         lsqr_res
%         relres(i)
%         lambda
        lambda = abs((etaeps-lsqr_res)/(relres(i)-lsqr_res))*lambda;
        lambdahis(i+1) = lambda; 
    end
    

    
    if strcmp(svdbasis,'x') 
        Xtemp = reshape(x,nim,nim);
        [UX,SX,VX] = svd(Xtemp);
        sx = diag(SX); 
        precX = kron(ones(nim,1),sx.^2);
        precX = precX + gamma;
        precX = precX.^((2-p)/2);
        gamma = gamma/2;
    end


    relerr(i) = norm(xtrue - X(:,i))/normxtrue;
    
    
    
end        
        
% relres = [beta1;relres];
% relerr = [norm(xtrue-x0)/normxtrue;relerr];