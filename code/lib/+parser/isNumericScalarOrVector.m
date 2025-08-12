function flag = isNumericScalarOrVector(x)
%function to validate input attributes used in inputParser: return true/false
%Check if input is a numberiv vector (1D, but could also be a scalar in degraded case)

flag = isnumeric(x) && isvector(x); 

end