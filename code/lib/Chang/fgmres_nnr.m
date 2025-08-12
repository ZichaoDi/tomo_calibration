function [xapprox,RelErr,RelRes,lambdahis,Zmatr,Vmatr,H] = fgmres_nnr(A, b, xex, x0,parameters)

m = parameters.maxIt;
reg = parameters.reg;
regmat = parameters.regmat;
svdbasis = parameters.svdbasis;
eta = parameters.eta;
nl = parameters.nl;
p = parameters.p;
thr = parameters.thr;
gamma = parameters.gamma;

N = length(b); n = sqrt(N);

xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);

nb = norm(b);
nx = norm(xex);


% B = reshape(b,n,n);
% [UB,~,VB] = svd(B);

% UX = UB; VX = VB; % in the case x0 = 0;
UX = eye(n); VX = eye(n);
SX = eye(n); 
%precX = SX(:);
sx = diag(SX); 

if strcmp(svdbasis,'x')
    precX = kron(ones(n,1),sx.^2);
    precX = precX + gamma;
    precX = precX.^((2-p)/4);
end

Vmatr = zeros(N,m+1);
Zmatr = zeros(N,m);
res = b - A*x0;
Vmatr(:,1)=res/norm(res);
nr = norm(res);

etaeps = eta*nl;

if strcmp(reg,'discrep')
    lambda = 1e8;
    % lambda = 1;
elseif isnumeric(reg)
    lambda = reg;
end
lambdahis = ones(m+1,1);

for k = 1:m
    
    
%%%%
  vtemp = Vmatr(:,k);
  vtemp = reshape(vtemp, n, n);

  
  if strcmp(svdbasis,'v')
      [UV, SV, VV] = svd(vtemp);
%       sv = diag(SV);

% %       diffsV = sv(1:end-1)-sv(2:end); 
% %       diffsVr = diffsV./diffsV(1);
% %       trunc = find(diffsVr<thr, 1);
% %       sv(trunc+1:end) = 0;

%       precX = kron(ones(n,1), sv);
%       precX = precX.^((2-p)/2);
%       vtemp = UV'*(vtemp*VV);
%       vtemp = precX.*vtemp(:);

      precX = SV.^(2-p/2);
      vtemp = precX(:);
      vtemp = reshape(vtemp, n, n); 
      vtemp = UV*(vtemp*VV');
  else
      if k>1
      vtemp = UX'*(vtemp*VX);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, n, n); 
      vtemp = UX*(vtemp*VX');
      end
  end
  Zmatr(:,k) = vtemp(:);
  
  
%%%%

    w = A*Zmatr(:,k);
    
    for i=1:k
        H(i,k)=w'*Vmatr(:,i);
        w=w-H(i,k)*Vmatr(:,i);
    end
    H(k+1,k)=norm(w);
    if H(k+1,k)<1e-10
        disp('Breakdown of GMRES')
        return
    else
        
        Vmatr(:,k+1)=w/H(k+1,k);
        d=[nr;zeros(k,1)];
        
        [UH,SH,VH] = svd(H(1:k+1,1:k));
        gmres_res = abs((UH(:,k+1)'*d))/nb;

        if strcmp(regmat,'I')
            c = UH'*d;

            if k == 1
                SH = SH(1,1);
            else
                SH = diag(SH);
            end

            Filt = SH.^2 + lambda;
            yj = VH*((SH.*c(1:k))./Filt);
            
        elseif strcmp(regmat,'R')
            [~,R] = qr(Zmatr(:,1:k),0);
            yj = [H(1:k+1,1:k);sqrt(lambda)*R]\[d;zeros(k,1)];
        end
        
        RelRes(k) = norm(H(1:k+1,1:k)*yj - d)/nb;
        
        if strcmp(reg,'discrep') || (isstruct(reg) && strcmp(reg.rule,'discrep'))
            lambda = abs((etaeps-gmres_res)/(RelRes(k)-gmres_res))*lambda;
            % lambdahis(k+1) = lambda; 
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
        
%         if strcmp(reg,'discrep')
%             beta = 1/lambda;
%             SH = [diag(SH); zeros(1,k)];
%             matrtemp = (beta*(SH*SH') + eye(k+1));
%             zetabeta = matrtemp\c;
%             wbeta = matrtemp\zetabeta; 
%             
%             f = (RelRes(k))^2;
%             f1 = 2/beta*zetabeta'*(wbeta - zetabeta); f1 = f1/(nr^2);
%             
%             % f1 = ...
%             beta=beta-(f-nl^2)/(f1); 
%             lambda = 1/beta;
%             
% %             lambda = abs((etaeps-gmres_res)/(RelRes(k)-gmres_res))*lambda;
%             lambdahis(k+1) = lambda; 
%       
%         end
        
    end
    xj = x0 + Zmatr(:,1:k)*yj;
    
    xapprox(:,k) = xj(:);
    RelErr(k) = norm(xj-xex)/nx;
    

    if strcmp(svdbasis,'x') 
        Xtemp = reshape(xj,n,n);
        [UX,SX,VX] = svd(Xtemp);
        sx = diag(SX); 
        precX = kron(ones(n,1),sx.^2);
        precX = precX + gamma;
        precX = precX.^((2-p)/4);
        gamma = gamma/2;
    end


end

% RelRes = [nr;RelRes];
% RelErr = [norm(xex-x0)/nx;RelErr];