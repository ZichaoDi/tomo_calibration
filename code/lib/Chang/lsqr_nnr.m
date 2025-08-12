function [xapprox,RelErr,RelRes,lambdahis,Vmatr,Bmatr] = lsqr_nnr(A, b, xex, x0,parameters)

m = parameters.maxIt;
reg = parameters.reg;
eta = parameters.eta;
nl = parameters.nl;
p = parameters.p;
UX = parameters.svd.U;
VX = parameters.svd.V;
SX = parameters.svd.S;
gamma = parameters.gamma;
% precX = SX(:);
% precX = precX.^((2-p)/2);

M = length(b); N = length(x0);
nim = sqrt(N);

sx = diag(SX); 
precX = kron(ones(nim,1),sx.^2);
precX = precX + gamma;
precX = precX.^((2-p)/4); % this is a vector

xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);

nb = norm(b);
nx = norm(xex);


Bmatr = zeros(m+1,m);
Umatr = zeros(M,m+1);
Vmatr = zeros(N,m+1);
Zmatr = Vmatr;
r = b - A_times_vec(A,x0);
beta = norm(r); nr = beta;
Umatr(:,1) = r/beta;


etaeps = eta*nl;

if strcmp(reg,'discrep') && parameters.cycle == 1
    lambda = 1;
elseif strcmp(reg,'discrep') && parameters.cycle > 1
    lambda = parameters.lambda1;
elseif strcmp(reg, 'opt')
    lambda = 0;
elseif isnumeric(reg)
    lambda = reg;
elseif isstruct(reg)
    % rule = reg.rule;
    lambda = reg.value;
    stagn = 0;
end

lambdahis = ones(m+1,1); % history of the lambda's
lambdahis(1) = lambda;

for k = 1:m
    
    w = Atransp_times_vec(A, Umatr(:,k));  
    
    % precond transpose
    w = reshape(w,nim,nim);
    w = UX'*(w*VX);
    w = precX.*w(:);
%     w = reshape(w,nim,nim);
%     w = UX*(w*VX');
%     w = w(:);
   
    if k>1
        w = w - beta*Vmatr(:,k-1);
    end
    
    for jj = 1:k-1
        w = w - (Vmatr(:,jj)'*w)*Vmatr(:,jj);
    end
    
    alpha = norm(w);
    Vmatr(:,k) = w/alpha;
    u = Vmatr(:,k);
    
    % precond
%     u = reshape(u,nim,nim);
%     u = UX'*(u*VX);
    u = precX.*u(:);
    u = reshape(u,nim,nim);
    u = UX*(u*VX');
    
    Zmatr(:,k) = u(:);
    
    u = A_times_vec(A,Zmatr(:,k));
    u = u - alpha*Umatr(:,k);
    
    % reorthogonalize 
    for jj = 1:k-1
        u = u - (Umatr(:,jj)'*u)*Umatr(:,jj);
    end
    
    beta = norm(u);
    Umatr(:,k+1) = u/beta;
    Bmatr(k,k) = alpha;
    Bmatr(k+1,k) = beta;

    if alpha < 1e-10 || beta < 1e-10
        disp('Breakdown of GKB')
        return
    else
        
        d=[nr;zeros(k,1)];
        
        [UB,SB,VB] = svd(Bmatr(1:k+1,1:k));
        c = UB'*d;
        
        lsqr_res = abs(c(k+1))/nb;


        if k == 1
            SB = SB(1,1);
        else
            SB = diag(SB);
        end
        
        if strcmp(reg,'opt')
            xexhat = reshape(xex, nim, nim);
            xexhat = UX'*(xexhat*VX);
            xexhat = xexhat(:)./precX;
            xhat = VB'*(Vmatr(:,1:k)'*xexhat);
            sB = SB;
            myfun = @(lambda) norm(((conj(sB)./(abs(sB).^2+lambda)).*c(1:end-1) - xhat));
            options.TolX = eps;
            lambda = fminbnd(myfun, 0, sB(1),options);

        end
            

        Filt = SB.^2 + lambda;
        yj = VB*((SB.*c(1:k))./Filt);
            
        
        RelRes(k) = norm(Bmatr(1:k+1,1:k)*yj - d)/nb;
        
        if strcmp(reg,'discrep') || (isstruct(reg) && strcmp(reg.rule,'discrep'))
            lambda = abs((etaeps-lsqr_res)/(RelRes(k)-lsqr_res))*lambda;
            % lambdahis(k+1) = lambda; 
            lambdahis(k+1) = lambda; 
        elseif (isstruct(reg) && strcmp(reg.rule,'discrepvar')) % && k>1
%             RegParamk = fzero(@(l)discrfcn(l, Mk, ZRksq, rhsk, eta*NoiseLevel), [0, 1e10]);
%             RegParamVect(k) = RegParamk;
            betaold = 1/lambda;
            SH = [diag(SH); zeros(1,k)];
            matrtemp = (betaold*(SH*SH') + eye(k+1));
            zetabeta = matrtemp\c;
            wbeta = matrtemp\zetabeta; 
            
            f = (RelRes(k))^2;
            f1 = 2/betaold*zetabeta'*(wbeta - zetabeta); f1 = f1/(nr^2);
            
            % f1 = ...
            beta=betaold-(f-nl^2)/(f1); 
            lambda = 1/beta;
            if lambda <= 1e-8
                lambda = 1/betaold;
                if lambda <=1e-14
                    stagn = stagn + 1;
                end
            end
            if stagn > 4
                break
            end
        end
        lambdahis(k+1) = lambda; 
    end
    
%     z = Vmatr(:,1:k)*yj;
    
%     % precond 
%     z = reshape(z,nim,nim);
%     z = UX'*(z*VX);
%     z = precX.*z(:);
%     z = reshape(z,nim,nim);
%     z = UX*(z*VX');
    
        
    xj = x0 + Zmatr(:,1:k)*yj;
    
    xapprox(:,k) = xj(:);
    RelErr(k) = norm(xj-xex)/nx;
    
%     if strcmp(reg,'opt')
%         lambdahis = lambdahis(2:end);
%     end
    
    % Xtemp = reshape(xj,nim,nim);
    
end

