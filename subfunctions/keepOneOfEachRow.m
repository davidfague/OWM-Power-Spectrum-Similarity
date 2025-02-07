function tblOut = keepOneOfEachRow(tblIn)
% KEEPONEOFEACHROW removes extra duplicates but keeps one occurrence of each row.
%
%   tblOut = unique(tblIn, 'rows', 'stable');
%   This preserves the first occurrence of each unique row, in the same order
%   as they appeared in tblIn.

    % Use unique to keep exactly one copy of each distinct row
    tblOut = unique(tblIn, 'rows', 'stable');
end
