function [x, s] = calcFeasibleXS(x0, A, b, tau)
%compute (x,s) > 0 that satisfies dual feasiblility (KKT stationary) condition s = A'*(Ax-b) + tau*1, by update x0 so that
%x = x0 + alpha*dx > 0, where dx = 1
%s = s0 + alpha*ds = A'*(Ax-b) + tau*1 = A'*(A*(x0+epsilon)-b) + tau*1 = (A'*(A*x0 - b) + tau*1) + alpha*A'*A*dx > 0
% we assume that dx, ds > 0, which is at least physically true for MFM
% Inputs: 
%       x0: Starting solution, it could be approximatley close to the solution
%       A:  An @A_operator that containts both AFun and AtFun operator

x = x0(:);
N = numel(x);

dx = ones(N,1);
s = A'*(A*x-b) + tau;
ds = A'*(A*dx); % typically close to all one vector  ,

alpha  = 0;
zAll = {[x, dx],  [s, ds]};
for i = 1:2
    z = zAll{i}(:,1);
    dz = zAll{i}(:,2);
    iz = (z<=eps);
    if nnz(iz)>0
        % find smallest epsilon_z, so that z + epsilon_z*dz >= smallPositive > 0,
        zNorm = norm(z,2);
        if zNorm < eps*N % z is all zero
            smallPositive = 1e-8;
        else
            smallPositive = zNorm/N*1e-8;
        end
        alpha_z = max((smallPositive - z) ./ dz);
        alpha = max(alpha, alpha_z);
    end
end

x = x + alpha * dx;
s = s + alpha * ds;

end

