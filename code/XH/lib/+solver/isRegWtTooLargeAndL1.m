function ret = isRegWtTooLargeAndL1(opts)
%solver.isRegWtTooLargeAndL1(opts), when regWt is too large, just need to return a solution of all 0
% Case I: 
% x_opt = argmin 0.5*||Ax-b||_2^2 + regWt*||x||_1
% x_opt = 0 if regWt >= max(abs(A'b)) derived from the condition that 0 belongs to its subgradient
% Case II:
% x_opt = argmin 0.5*||Ax-b||_2^2 + regWt*||x||_1 = argmin 0.5*||Ax-b||_2^2 + regWt*1'*x
%       s.t. x >= 0
% x_opt = 0 if regWt >= max(A'b) derived from the KKT condition



%% L1 norm and soft threshold function
L1Fun = @(x) sum(abs(x(:)));                % L1 norm,
%softFun = @(x,T) sign(x).*max(abs(x)-T, 0); % soft threshold for real. If complex, b = max(abs(x) - T, 0); b = b./(b+T) .* x;

Atb = opts.AtFun(opts.b);
ret = (opts.regWt >= max(abs(Atb(:)))) ...
    && (isequal(opts.regFun, L1Fun) ||  opts.regFun(Atb)==L1Fun(Atb));
end


% ret = (opts.regWt >= max(abs(Atb(:)))) ...
%     && (isequal(opts.regFun, L1Fun) ||  opts.regFun(Atb)==L1Fun(Atb)) ...
%     && (isequal(opts.denoiseFun, softFun) || isequal(opts.denoiseFun(Atb, opts.regWt), softFun(Atb, opts.regWt)));
% 
