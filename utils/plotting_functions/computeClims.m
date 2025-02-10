function clims = computeClims(data)
% computeClims Computes color limits for imagesc based on the data's statistics.
%
%   clims = computeClims(data) returns a two-element vector [min, max] where:
%       min = mean(data) - 2*std(data)
%       max = mean(data) + 2*std(data)
%
%   This is useful for setting the color limits in an imagesc plot, e.g.:
%
%       data = randn(100,100);
%       imagesc(data, computeClims(data));
%       colorbar;
%
%   The function treats the data as a single vector by using data(:) to
%   compute the overall mean and standard deviation.
%
%   Inputs:
%       data - A numeric array (e.g., matrix) that contains the data to be displayed.
%
%   Outputs:
%       clims - A two-element vector [min, max] for use as the color limits in imagesc.
%

    % Compute the mean and standard deviation over all elements of the data
    mu = mean(data(:));
    sigma = std(data(:));
    
    % Define the color limits as 2 standard deviations around the mean
    clims = [mu - 2*sigma, mu + 2*sigma];
end
