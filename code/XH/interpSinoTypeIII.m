function SInterp = interpSinoTypeIII(S, drift, interpMethod)
%interpSinoTypeIII polynomial interpolation of sinogram
% Input:   
%       S:      2D original sinogram of size NTau*NTheta, without drift
%       drift:  Type III drift. NTheta*NTau, assume small mostly [-0.5, 0.5], TODO: add check if it is too large
%       interpMethod:    interplolation options
%               'linear':   
%               'central_diff': use S(i,j-1) and S(i,j+1) to interpolate result
%               'quadratic' etc.: not inploated , 
% Output:   
%       SInterp: 2D interpolated sinogram of size NTau*NTheta, which is an approximation of S drifted by drift
% See also:
%       calcDriftTypeIII(); interpSino();


%%
if nargin < 1; testMe(); end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'linear';
end


[NTheta, NTau] = size(S);

%%
%SInterp = zeros(size(S), 'like', S);

% 1st order difference: forward/central difference:
SDiffFwd = zeros(size(S),'like',S); n = 1:(NTau-1); SDiffFwd(:,n) = S(:,n+1) - S(:,n); % SDiffFwd = [S(:, 2:end) - S(:, 1:end-1), zeros(NTheta, 1)]; 

% compute integer and fractional part, which is used for linear/quadratic opitions
d = floor(drift);  alpha = drift - d;  % 0 <= alpha < 1

switch lower(interpMethod)
    case 'linear'
        % At angle theta, To recover beam i of SInterp, we use two cloest beams j and j+1 (of S) to linear interpolate, where j = Tau2(theta, i).
        % SInterp(theta, i) =  S(theta, j) + alpha(theta, i)*(S(theta,j+1) - S(theta,j)); % 0 <= alpha(i) < 1 is the drifted distance to beam j (in S) from beam i (in SInterp)        
        % 2D beams: SInterp(theta, tau) uses S(theta, tau2), S(theta, tau2+1) to interpolate, where tau2=Tau2(theta, tau)
        [Tau, Theta] = meshgrid(1:NTau, 1:NTheta);
        Tau2 = Tau + d;  Tau2(Tau2<1) = 1;  Tau2(Tau2>NTau-1) = NTau-1; 
        
        % vectorized version
        idx = sub2ind(size(S), Theta, Tau2);        
        SInterp = S(idx) + SDiffFwd(idx) .* alpha;
        
        % loop version. Only keep for reference and verification, same effect but slower than above vectorized version
%         for theta = 1:NTheta
%             for tau = 1:NTau
%                 tau2 = Tau2(theta, tau);
%                 SInterp(theta, tau) = S(theta, tau2) + SDiffFwd(theta, tau2)*alpha(theta, tau);                
%             end
%         end
    %case 'quadratic'  
    case 'central_diff'
        % test interpolation and reconstruction, assume drift is -1 to 1, and use pixels left and right for central difference
        SDiffCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SDiffCtr(:,n) = (S(:,n+1) - S(:,n-1))/2; % SDiffCtr = [zeros(NTheta, 1), (S(:, 3:end) - S(:, 1:end-2))/2, zeros(NTheta, 1)];
        % SMeanCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SMeanCtr(:,n) = (S(:,n+1) + S(:,n-1))/2; 
        SMeanCtr = S; 
        SInterp = SDiffCtr.*drift + SMeanCtr;
    otherwise
        error('Unknown interpolation method: %s.', interpMethod);        
end



end


function testMe(SAll, driftGT)
%% 2020.04.03_verfied that the interpSinoTypeIII is good at least for 50x50 Phantom with sinogram 45x74, drift std_=0.5
S = SAll(:,:,1);
SInterpGT = SAll(:,:,2);
SInterp = interpSinoTypeIII(S, driftGT);
fprintf('The psnr(S, SInterpGT) %.2f dB\n', difference(S, SInterpGT)); %23.01 dB
fprintf('The psnr(SInterp, SInterpGT) %.2f dB\n', difference(SInterp, SInterpGT)); %29.39 dB

figure, multAxes(@imagesc, {S-SInterpGT, SInterp-SInterpGT}); linkAxesXYZLimColorView(); multAxes(@colorbar);multAxes(@title, {'S-SInterpGT', 'SInterp-SInterpGT'})
% error analysis, why some error are large, and some are small
figure, multAxes(@imagesc, {S, SInterpGT}); linkAxesXYZLimColorView(); multAxes(@colorbar);multAxes(@title, {'S', 'SInterpGT'})
[Grad_tau, Grad_theta] = imGrad(S);
figure, multAxes(@imagesc, {abs(Grad_theta), abs(Grad_tau)}); linkAxesXYZLimColorView(); multAxes(@colorbar);multAxes(@title, {'abs(Grad_theta)', 'abs(Grad_tau)'})

mask = (S>eps);
SS = SInterpGT; SS = SS + mean(SS(mask)); SS(SS<eps) = inf;
figure, multAxes(@imagesc, {(S-SInterpGT)./SS, (SInterp-SInterpGT)./SS}); linkAxesXYZLimColorView(); multAxes(@colorbar);multAxes(@title, {'(S-SInterpGT)./SS', '(SInterp-SInterpGT)./SS'})
driftMap = abs(driftGT); driftMap(~mask) = 0;
figure, multAxes(@imagesc, {driftMap, driftMap}); linkAxesXYZLimColorView(); multAxes(@colorbar);multAxes(@title, {'driftMap', 'driftMap'})



end
