function MInvFun = getMInvFun(preconditioner, d, rho, diagATA, Vk)
%getMInvFun() M for (DInv + A'*A)
%Usage
% MInvFun = solver.getMInvFun([]); or MInvFun = solver.getMInvFun('');  % no predonctioner, so MInvFun(x) = x;
% MInvFun = solver.getMInvFun(MInvFun); %MInvFun is already a function handle where MInvFun(x) compute M\x.
% MInvFun = solver.getMInvFun('zeroATA', dInv); % A'*A is approximate by zeros(N,N) so that M = DInv;
% MInvFun = solver.getMInvFun('identityATA', dInv, rho); % rho is a scalar so that is a A'*A is approximated by rho*eye(N,N);
% MInvFun = solver.getMInvFun('diagATA', dInv, [], diagATA); % diagATA is a vector so that A'*A is approximate by diag(diagATA);
% MInvFun = solver.getMInvFun('lowRankATA', dInv, [], [], v); % v is a vector so that A'*A is approximated by v*v
% MInvFun = solver.getMInvFun('lowRankATA', dInv, [], [], Vk); % Vk is N*k matrix so that A'*A is approximated by Vk*Vk'

%
% system (DInv + ATA)*dx = c
% can be solved using pcg with MInvFun as:
% pcg(H, c, tolCG, maxIterCG, MInvFun, [], x0)
%
% M is the approximate of A
% M = P^T*P is what we used in the paper AlgorithmReconstruct3D, we only provide MInvFun() the function handle MInvFun(x) = M\x;
% M = C^T*C is used in ``Numerical Optimization'', Jorge Norcedal and Stephen Wright, 2nd Edition, Page 113, where C = P
% --------------- In Matlab pcg()--------------
% M = M1*M2 is used matlab function pcg(),  where  M1 = P^T; M2 = P. NOTE: M=M1*M2 also used in http://netlib.org/linalg/html_templates/node54.html
% Note in matlab if the input argument is a matrix, then it is M; if it is a function handle MInvFun(), then MInvFun(r) returns M\r, it is actually function of M^(-1)
%
% Inputs:
%   preconditioner: preconditioner for PCG. It could be three formats function handle, empty matrix/strings, or strings
%       1. function handle:  preconditioner is already function handle of MInvFun(), just return itself
%       2. empty matrix [] or empty strings '':  This means we don't use preconditioner or preconditioner is an identy matrix, MInvFun(x) = x;
%       3. following strings:
%           3.1 'zeroATA':      approximate A'*A as a all zero matrix.
%           3.2 'identityATA':  approximate A'*A as a scaled identity matirx rho*ATA
%           3.3 'diagATA':      approximate A'*A as diagonal matrix
%           3.4 'lowRankATA':   approximate A'*A as rank 1 matrix lambda*v*v' use its largest eigenvalue lambda with corresponding eigenvector v%
%  d:       a vector d = diag(D) = 1./diag(DInv)
%  rho:     used for 'identityATA' option, where [M,N]=size(A), here A'*A is approximate by rho*eye(N, N);
%  diagATA: used for 'diagATA' option, where A'*A is approximate by diag(diagATA), where diagATA should be diag(A'*A).
%  Vk:      used for 'lowRankATA'option, where A'*A is approximate by V*V', where V is N*K matrix and K << N;
%                               V is just the first K orthnormal eigen vector scaled by square root of eigen vector.
%                               V = V0(:, 1:K)*sqrt(Lambda(1:K, 1:K)), where V0 is the orthnormal eigen vector matrix, Lambda is
%                               diagonal matrix whose diagonal elements are eigen values sorted descendly, V0*Lambda*V0' = A'*A
%                               from [V0, Lambda] = evd(A'*A);
%                               For example when K = 1, V is vector v, and A'*A is approximated by rank 1 matrix
%                               v*v'=lambda*v0*v0', where lambda is the largest eigenvalue and v0 is the corresponding eigenvector
%
% NOTE: In MFM, (AAT)(x) can be much faster than A(AT(x)), but (ATA)(x) can only computed simlar speed as (AT


% 1. function handle
if isa(preconditioner, 'function_handle')
    MInvFun = preconditioner;
    return
end

% 2. empty matrix
if isempty(preconditioner)
    MInvFun = @(x) x;
    return
end

% 3. strings
assert(isstr(preconditioner), 'The preconditioner variable must be the following three format: function handle, empty matrix/strings, or strings with meanings');
switch lower(preconditioner)
    case lower('zeroATA')
        MInvFun = @(x) d.*x; % same as x./dInv;
    case lower('identityATA')
        dInv = 1./d;
        HInvApprox = 1./ (dInv + rho);
        MInvFun = @(x) HInvApprox.*x;
    case lower('diagATA')
        dInv = 1./d;
        HInvApprox = 1./ (dInv + diagATA);
        MInvFun = @(x) HInvApprox.*x;
    case lower('lowRankATA')
        [~, k] = size(Vk);
        if k == 1
            % special rank 1 case : it can be computed use the same formula as general case below, but I like write this out for vector, no inverse, or \ just dividion.
            v = Vk;
            scale = 1/(1 + v'*(d.*v));
            MInvFun = @(x) d.*(x - v*(scale*(v'*(d.*x)))); %MInvFun = @(x) d.*(x- v*(v'*(d.*x)/(1+v'*(d.*v))));
        else
            % general rank k case:
            SCALE = inv(eye(k, k) + Vk'*bsxfun(@times, Vk, d)); % (I+Vk'*D*Vk)^{-1}
            MInvFun = @(x) d.*(x - Vk*(SCALE*(Vk'*(d.*x))));
        end
    otherwise
        error('unkown precondtioner named: %s', preconditioner);
end
end
