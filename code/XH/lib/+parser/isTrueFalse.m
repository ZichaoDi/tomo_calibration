function flag = isTrueFalse(x)
%function to validate input attributes used in inputParser: return true/false
%Check if input is one of true/false/1/0

flag = isscalar(x) && ismember(x, [true, false, 1, 0]); % true/false/1/0

end