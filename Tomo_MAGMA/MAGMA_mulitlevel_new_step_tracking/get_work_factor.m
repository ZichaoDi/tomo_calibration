% get_work_factor.m - Determines the scaling factor for a given MAGMA level.
function factor = get_work_factor(level)
    % The work of one iteration at a given level is scaled relative to a fine-level iteration.
    % Level 1 (Fine): Factor = 1
    % Level 2 (Coarse): Factor = 1/4
    % Level 3 (Coarser): Factor = 1/16
    % Level L > 1: Factor = 1 / (4^(L-1))
    if level == 1
        factor = 1.0;
    else
        factor = 1 / (4^(level - 1));
    end
end