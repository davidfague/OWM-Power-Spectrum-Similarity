function frequencies = get_frequencies(num_frequencies, range, spacing)
%   frequencies = get_frequencies(num_frequencies, range, spacing) returns a vector of 
%   num_frequencies equally spaced between range(1) and range(2).
%
%   - num_frequencies: number of points to generate.
%   - range: a two-element vector [min, max] specifying the range.
%   - spacing: a string that can be 'linear' or 'log' to specify the type
%     of spacing. Default is 'linear' if not provided.
%
%   Example usage:
%      x = get_frequencies(50, [1, 100], 'log');    % Logarithmically spaced points
%      y = get_frequencies(100, [0, 10], 'linear');   % Linearly spaced points

    if nargin < 3
        spacing = 'linear';
    end

    % Validate the range vector
    if numel(range) ~= 2 || range(1) >= range(2)
        error('Range must be a two-element vector [min, max] with min < max.');
    end

    switch lower(spacing)
        case 'linear'
            frequencies = linspace(range(1), range(2), num_frequencies);
        case {'log', 'logarithmic'}
            % For logarithmic spacing, ensure the starting value is positive
            if range(1) <= 0
                error('For logarithmic spacing, the starting point must be > 0.');
            end
            frequencies = logspace(log10(range(1)), log10(range(2)), num_frequencies);
        otherwise
            error('Unknown spacing type. Use ''linear'' or ''log''.');
    end
end