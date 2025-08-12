function [AFun, AtFun, B] = getAFunAtFunB(Sm, S0_rec, maxDrift)
% Twist-linear:  0.5*||A*D(:) - B(:)||.^2 +  regWt*regFun(D) 
% Twist-nonlinear:  0.5*||A(D(:)) - B(:)||.^2 +  regWt*regFun(D) 
% Twist-linear: A is diagonal matrix for tomography distortion. Each scan output does NOT depend on other scans locations
% Twist-nonlinear: A() is a vector of scalar functions
%   S:  the Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
% maxDrift is interger limit 


if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

if ndims(S0_rec) == 2
    [NTheta, NTau] = size(S0_rec);
    S_sz = [NTheta, NTau];

    SDiffCtr = zeros(size(S0_rec),'like',S0_rec); 
    n = (maxDrift+1):(NTau-maxDrift); 
    SDiffCtr(:,n) = (S0_rec(:,n+maxDrift) - S0_rec(:,n-maxDrift))/2; 
    
    AFun = @(D) reshape(SDiffCtr(:).*D(:), S_sz); % AFun  = @(D) reshape(A*D(:), S_sz); % where D is 2D drift matrix
    AtFun = AFun;
    
    SMeanCtr = zeros(size(S0_rec),'like',S0_rec); 
    SMeanCtr(:,n) = (S0_rec(:,n+maxDrift) + S0_rec(:,n-maxDrift))/2;     
    B = Sm - SMeanCtr; 
    
    % B = Sm - S0_rec;  % different version, use all S+, S- and S0, instead of use mean of S+/S-, use S0, result similar
    
else
    %%
    assert(ndims(S0_rec)==3, 'dimension must be 2 or 3')
    %S0_rec is actually S_resAll, super-res version
    S_resAll = S0_rec;
    %%
    [NTheta, NTau, N_res] = size(S_resAll); N_res  = (N_res + 1)/2
    
end
