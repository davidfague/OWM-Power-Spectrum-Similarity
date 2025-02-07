function out = removeDirectionPrefix(strs)
% removeDirectionPrefix removes a leading "l", "r", "left", or "right"
% (case insensitive) and any spaces immediately following that token from
% each string in the input.
%
%   out = removeDirectionPrefix(strs)
%
% Input:
%   strs - A column (or vector) of strings. This can be a string array or
%          a cell array of character vectors.
%
% Output:
%   out  - A column (or vector) of strings with the specified prefix removed.
%
% Example:
%   strs = ["Left apple", "right banana", "l orange", "R mango"];
%   out = removeDirectionPrefix(strs)
%   % out is:
%   %    "apple"
%   %    "banana"
%   %    "orange"
%   %    "mango"

    % Define a regular expression pattern:
    %   ^           : match the start of the string
    %   (l(?:eft)?  : match 'l' optionally followed by 'eft'
    %    |r(?:ight)?): or match 'r' optionally followed by 'ight'
    %   \s*         : match any number of whitespace characters that follow
    pattern = '^(l(?:eft)?|r(?:ight)?)\s*';

    % Use regexprep to replace the matching pattern with an empty string.
    % The 'ignorecase' option makes the matching case-insensitive.
    out = regexprep(strs, pattern, '', 'ignorecase');
end