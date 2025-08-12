function sampledHist = updateSampledHist(sampledHist, x, cpuTimes, k, sampleInterval)
% sampledHist is a structure array that keeps history of x, and corresponding k, cpuTime at every sampleInterval
% Two cases of sampleInterval
% Case I:   sampleInterval is a scalar, sampleInterval is the time interval in seconds, a iteration so far takes interger number 
%           of sampleInterval time will be saved, such as sampleInterval = 3600 keeps a iteration every 3600 seconds.
% Case II:  sampleInterval is a vector(1D array), sampleInterval is the iteration numbers, result of those iteration numbers will be saved
%           Here we could do non-uniform sampled, such as let sampleInterval = [100, 200, 1000] to keep 100/200/1000 iterations
% Case III: sampleInterval is a cell array (2x1, 1x2 or more), then each cell element specifies keep rules either in Case I or II, 
%           and keep current history if any of this cell element.
%           sampleInterval = {3600, [100, 200, 500, 1000]} will both keep iteration every 3600s and iteration#[100, 200, 500, 1000
% NOTE: k-1 is the iteraiton number. k=2 means the first iteration. k=1 means the initial (0th) iteration.


%% Parse and verify Inputs
if ~exist('sampleInterval', 'var') || isempty(sampleInterval)
    sampleInterval = 3600; % default 1 hour = 3600 seconds
end

assert(k>=1, 'Only save sampled history for k>=1. Ignore k=0 which is the initial!');
%%    

if iscell(sampleInterval)
    % Case III
    isKeep = false;
    for n = 1:numel(sampleInterval)
        isKeep = isKeep || calcIsKeep(sampleInterval{n}, cpuTimes, k);
    end
else
    % Case I or II
    isKeep = calcIsKeep(sampleInterval, cpuTimes, k);
end

if isKeep
    % k+1 is the first iteration that equal or more than integer muliplies of histTimeInterval.
    % Edge case cpuTimes(k+1) = cpuTimes(k) = Interger*histTimeInterval is also handled here, though in practice cpuTime
    % should always be xxx.xxxx seconds and cpuTimes(k+1) = cpuTimes(k) may never happen
    curStep = struct('x', x, 'k', k, 'cpuTime', cpuTimes(k+1));
    if isempty(sampledHist) || (isstruct(sampledHist) && numel(sampledHist)==1  && all(cellfun(@(x) isempty(sampledHist.(x)), fieldnames(sampledHist))))
        % keep the first history, when sampledHist == [] or sampledHist =  struct('x', [], 'k', [], 'cpuTime', []);
        sampledHist = curStep;
    else
        sampledHist(end+1) = curStep;
    end    
end

end

function isKeep = calcIsKeep(sampleInterval, cpuTimes, k)
if isscalar(sampleInterval)
    % Case I:
    isKeep = floor(cpuTimes(k+1)/sampleInterval) >= floor(cpuTimes(k)/sampleInterval) + 1;
else
    % Case II: 
    isKeep = any(k-1 == sampleInterval); % iteration number is k-1
end
end
