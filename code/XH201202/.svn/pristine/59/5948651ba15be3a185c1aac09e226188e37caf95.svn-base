function [obj, grad] = sfun_drift(d,w, S_measured, WSz)
    global L_d LNormalizer
% 2020.09.22: for a given w, d, compute the gradient of g(w, d) with respect to [w, d]
% assume using +1, -1 grid
% L_d = L(d) the forward model matrix given drift
% s_simulated = s(w, d) in paper, S_simulated is the 2D versoin
% s_measured = s^m in paper, S_measured is the 2D version
[NTheta, NTau] = size(S_measured);
SSz = [NTheta, NTau]; 

N = prod(WSz);
L_d = XTM_Tensor_XH(WSz, NTheta, NTau, reshape(d, SSz)); %WGT
L_d = L_d/LNormalizer;
S_simulated = reshape(L_d * w, SSz);

r = S_simulated(:) - S_measured(:);

% SDiffCtr(i) = a_i = dsi/ddi in paper
SDiffCtr = zeros(size(S_simulated),'like', S_simulated); n = 2:(NTau-1); SDiffCtr(:,n) = (S_simulated(:,n+1) - S_simulated(:,n-1))/2; 
grad = r.* SDiffCtr(:); 
obj = r'*r/2;
