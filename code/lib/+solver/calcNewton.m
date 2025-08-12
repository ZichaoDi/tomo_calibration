function [dx, itersCGNo, cgRelRes] = calcNewton(A, d, c, tolCG, maxIterCG, diagATA, isPreconditioned, isADAT, dx0, dy0, MInvFun)
%2018.02.02_depreciate, pcg() with solver.calcPrecondioner() solves the purpose well
%
% use MInvFun = solver.calcPrecondioner() to compute the precondioner, and call matlab
%Calculate Newton Step dx from equation:  H*dx = c
% where H = ATA + DInv
% where DInv = diag(dInv) is a diagonal matrix
%CALCNT Summary of this function goes here
%   Detailed explanation goes here
% NOTE: In MFM, (AAT)(x) can be much faster than A(AT(x)), but (ATA)(x) can only computed simlar speed as (AT
% dx0/dy0:  Intial value used for cg;  showed use dx/dy from last step performance very close to use all zero as initial.
%
% preconditioner M is the approximate of A
% M = P^T*P is what we used in the paper AlgorithmReconstruct3D, we only provide MInvFun() the function handle MInvFun(x) = M\x;
% M = C^T*C is used in ``Numerical Optimization'', Jorge Norcedal and Stephen Wright, 2nd Edition, Page 113, where C = P
% --------------- In Matlab pcg()--------------
% M = M1*M2 is used matlab function pcg(),  where  M1 = P^T; M2 = P. NOTE: M=M1*M2 also used in http://netlib.org/linalg/html_templates/node54.html
% Note in matlab if the input argument is a matrix, then it is M; if it is a function handle MInvFun(), then MInvFun(r) returns M\r, it is actually function of M^(-1)
% So we should input call matlab pcg with MInvFun as:  pcg(H, c, tolCG, maxIterCG, MInvFun, [], x0)


N = numel(c);
if ~exist('dx0', 'var') || isempty(dx0)
    dx0 = zeros(N,1);
end

if ~exist('dy0', 'var')
    dy0 = [];
end

dInv = 1./d;
HFun = @(dx) A'*(A*dx) + dInv.*dx; % function handle of (Dinv + A'A)*dx
%isPreconditioned = 1;
if exist('isADAT', 'var') && isADAT
    % woodbury inv(A'A+Dinv) = D-DA'inv(I+ADA')AD,
    % dx = (D-DA'inv(I+ADA')AD)*c = Dc - DA'*inv(I+ADA')*ADc = Dc - DA'*dy
    % solve (I + A*D*A')*dy = ADc.  NOTE: the diag(I+ADA') is tricky to compute, so no diagonal preconditionning here
    ADATIFun = @(dy) dy + A*(d .* (A'*dy));
    c_y = A*(d.*c);
    [dy, itersCGNo, cgRelRes] = cg(ADATIFun, c_y, dy0, tolCG, maxIterCG); % cg(A, b, x0, tol, maxIter), use previous dx as intial guess
    dx = d.*(c - A'*dy);
else
    % todo: change isPreconditioned to preconditioner with choice {'diagATA', 'identityATA', 'rank1ATA', MInvFun}; where MInvFun is the
    % function handle user specified, while the other three corresponding to three difference preconditioners
    if isPreconditioned
        % Solve (A'A+Dinv)*dx = c with preconditioner P = P^T; M = P^(-T)*P^(-1) = P^(-2);
        % P(A'A+Dinv)P*(Pinv*dx) = P*c
        % page 45:  1994_Shewchuk_An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
        pInv = sqrt(diagATA + dInv); % pInv^2 is the diag of hessianJ = A'A + Dinv
        p = 1./pInv; % p is the diagonal of the preconditional matrix P, where P^2 =  1./diag(hessianJ)
        % function of P(A'A+Dinv)P*dz, where dz = PInv*dx,  vector out for vector in: dz
        hPHP = @(dz) p.*(A'*(A*(p.*dz)) + dInv.*(p.*dz));
        %[dz, itersCGNo, cgRelRes] = cg(hPHP, p.*c, pInv.*dx0, tolCG, maxIterCG); % cg(A, b, x0, tol, maxIter), use previous dx as intial guess
        %dx = p.*dz;

%         [dz, flag, cgRelRes, itersCGNo] = pcg(hPHP, p.*c, tolCG, maxIterCG, [], [], pInv.*dx0); %%pcg(H,c,tol,maxit,M1,M2,x0); %matlab's own function, same result as our cg function
%         dx = p.*dz;
        p2 = 1./ (diagATA + dInv);
        MInvFunPCG =  @(x) p2.*x;
        %MInvFunPCG = @(x) MInvFun(x, d); % mainMFMTest with small case for #iteration in PCG: rankOne-PCG 28,  diagonal-PCG 99, no-conditioner CG 216;
        [dx, flag, cgRelRes, itersCGNo] = pcg(HFun, c, tolCG, maxIterCG, MInvFunPCG, [], dx0); % %pcg(H,c,tol,maxit,M1,M2,x0); %matlab's own function
        %difference(dx, dxNew);
        numel(p2);
    else
        % Directly Solve (A'A+Dinv)*dx = c.  NOT good! should never use this
        % Xiang's cg vs matlab pcg, cg converges much faster and results better accuracy due to the additional correction step
        % of r = b - A*x
        % cg:  Total iterations:  216, psnrSparse  =7.3182e+01dB
        % pcg: Total iterations:  715, psnrSparse =7.3075e+01dB
        %[dx, itersCGNo, cgRelRes] = cg(HFun, c, dx0, tolCG, maxIterCG); % cg(A, b, x0, tol, maxIter), use previous dx as intial guess
        [dx, flag, cgRelRes, itersCGNo] = pcg(HFun, c, tolCG, maxIterCG, [], [], dx0); % %pcg(H,c,tol,maxit,M1,M2,x0); %matlab's own function
    end
end

end
