function limits_to_use = get_clim_based_on_2z(data)
        flattened_data = data(:);
        % Calculate the mean and standard deviation of the flattened data
        data_mean = mean(flattened_data(:));
        data_std = std(flattened_data(:));
        % Set the color limits to be 2 standard deviations around the mean
        limits_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
end