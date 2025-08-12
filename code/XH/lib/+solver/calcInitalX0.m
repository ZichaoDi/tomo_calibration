function x0 = calcInitalX0(x0, AFun, AtFun, b, Atb)
% input x0 could be three cases:
%       (1) a string indicating the initalization choice, then return 
%       (2) [], no need to process, return []
%       (3) a numberic array, then no need to process, return x0, 

% Atb is AtFun(b), it is redundant but still passed here to save computation if it is avaliable from caller


if ~isempty(x0) && ischar(x0) 
    % case I: x0 is a string with choice of ininilazation    
    if ~exist('Atb', 'var') || isempty(Atb)
        Atb = AtFun(b);
    end
    xSz = size(Atb);
    switch lower(x0)
        case {'zero', 'zeros'}
            x0 = zeros(xSz, 'like', b);
        case {'random', 'rand'}
            x0 = rand(xSz, 'like', b);
        case lower('Atb')  % no! this is not good
            x0 = calcLsqInit(AFun, AtFun, b, Atb);  % typically is better than all zeros
        otherwise
            error('unknown tag ''%s'' to set x0 ', x0);
    end    
else
    % case II & III: x0 is [] or a numberic array
    assert(isempty(x0) || isnumeric(x0), 'Wrong x0 format! It need to be either [], or a string, or a numberic array!');    
end
end


function x0 = calcLsqInit(AFun, AtFun, b, Atb)    
%Simple initialization: x0 = scale* A'*b, 
% where scale is find from least-squares fit to minimize ||scale*A*A'*b - b||_2
% we could derive scale as
%      scale = b'*b2 / (b2'*b2) 
% where b2 = A*A'*b
% Atb is AtFun(b), it is redundant but still passed here to save computation if it is avaliable from caller

if ~exist('Atb', 'var') || isempty(Atb)
    Atb = AtFun(b);
end

b2 = AFun(Atb);
scale = (b(:)'*b2(:)) / (b2(:)'*b2(:));  % sum(AtFun(ones(size(b)))) = sum(A(:)) = sum(A'(:))
x0 = scale*Atb;

end

function x0 = calcRealLsqInit(AFun, AtFun, b, Atb)
% argmin ||Ax-b||^2 + epsilon*||x||^2 will get x0 = (A'A+epsilon*I)^(-1) A'b
% this could be solved for MFM after decompose A = F^H * O^H * O * F, which reduce to solve (O^H * O + epsilon*I)^(-1)
% tried to directly solve for small A
% always get half negative x for small epsilon; 
% only get all small positive numbers for large epsilon but if epsilon is large, then ||Ax-b||^2 is large which makes no sense
end


function x0 = calcLsqInitMeanMatchOld(AFun, AtFun, b, Atb)    
% to discard after verify the output is the same as calcLsqInit() (only verified in small cases not in real example), as the new function calcLsqInit() is much more clearer
%Simple initialization: x0 = scale* A'*b, 
% where scale is to make a least-squares fit to the mean intensity:  scale = mean(b)/mean(A*x0)
% change to 
% b(:) = A*f(:);
% x0 = A'*b(:) * normalize factor
% Atb is AtFun(b), it is redundant but still passed here to save computation if it is avaliable from caller

if ~exist('Atb', 'var') || isempty(Atb)
    Atb = AtFun(b);
end

f_bAllOnes = AtFun(ones(size(b)));
sumA = sum(f_bAllOnes(:)); % i.e., the summation of all elements in A matrix where b(:) = A*f(:);s
scale = (sum(b(:))*numel(Atb)) / (sum(Atb(:))*sumA);  % sum(AtFun(ones(size(b)))) = sum(A(:)) = sum(A'(:))
x0 = scale*Atb;

end
