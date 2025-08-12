function P = genDriftMatrix(drift, NTheta, NTau)
%Generate the drift matrix that maps undrifted model L0 to drifted model L
% Input: 
%   drift: the drift amount for each beamline (assume invariant with rotation angle),  NTau*1 vector
% Output: 
%   P:  NTheta*NTau by NTheta*NTau sparse matrix that L = P*L0, where L0 is the model matrix without drift, L is the model
%   matrix with drift
% We currently uses linear interplolation
% See also:
%       interpSino()
%       

d = floor(drift);
alpha = drift - d;

%P = sparse(NTheta*NTau, NTheta*NTau); % P = sparse(r, c, v, NTheta*NTau, NTheta*NTau); 
% row indices, column indices, value vector to specify P, there are at most two non-zeros in each column.
r = ones(NTheta*NTau*2, 1); 
c = ones(NTheta*NTau*2, 1); 
v = zeros(NTheta*NTau*2, 1); 
idxStart = 1; % for ith iteration, the ith block matrix is defineded by starting index of r(idx),c(idx),v(idx), where idx = idxStart +(1:NTheta)


for i = 1:NTau
    j = i+d(i); 
    
    %S2(:, i) = (1-alpha(i)) * S(:, j)  + alpha(i)*S(:, j+1);
    if (j >= 1) && (j <= NTau)
        %S2(:, i) = (1-alpha(i)) * S(:, j);
        idx = idxStart + (1:NTheta);
        r(idx) = (1:NTheta) + (i-1)*NTheta;
        c(idx) = (1:NTheta) + (j-1)*NTheta;
        v(idx) = 1 - alpha(i);
        idxStart = idxStart + NTheta;
        
        if (alpha(i) > eps) 
            if ((j+1) <= NTau)  % j >= 1 already, no need to check j+1 >= 1
                %S2(:, i) = S2(:, i) + alpha(i)*S(:, j+1);
                idx = idxStart + (1:NTheta);
                r(idx) = (1:NTheta) + (i-1)*NTheta;
                c(idx) = (1:NTheta) + j*NTheta;
                v(idx) = alpha(i);
                idxStart = idxStart + NTheta;                
            else                        
                error('Out of boundary! You should not reach here: %d, alpha(i)=%.2f.\n', d(i), alpha(i));
            end
        end        
        %NOTE: if alpha(i) < eps, then S(:, j+1) is not used, and we no need to check whether j+1 is within [1, NTau].
    else
        error('Out of boundary! You should not reach here: %d, alpha(i)=%.2f.\n', d(i), alpha(i));
    end
end


idxUsed = 1:idxStart-1;
P = sparse(r(idxUsed), c(idxUsed), v(idxUsed), NTheta*NTau, NTheta*NTau);  

end

