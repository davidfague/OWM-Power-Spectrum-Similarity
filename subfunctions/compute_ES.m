function [similarity_matrix, window1_mean_PS_vectors, window2_mean_PS_vectors] = compute_ES(Zpower, window_start_times, window_end_times, window_IDs1, window_IDs2, ES_freq_band, save_mean_PS_vectors)
    % Zpower is frequency X nTimes
    % window_start_times and window_end_times are same length
    % window_IDs select within start/end time arrays.
    % ES_freq_band ex. 1:40
    % for every window, get the average Zpower in the window
    % output matrix should be Nwindow_IDs1 x Nwindow_IDs2
    % mean_PS_vectors are Nwindow_IDs x Nfrequencies
    similarity_matrix = zeros(length(window_IDs1), length(window_IDs2));

    % Only allocate memory for mean PS vectors if needed
    if save_mean_PS_vectors
        window1_mean_PS_vectors = zeros(length(window_IDs1), length(ES_freq_band));
        window2_mean_PS_vectors = zeros(length(window_IDs2), length(ES_freq_band));
    else
        window1_mean_PS_vectors = [];
        window2_mean_PS_vectors = [];
    end

    for i = 1:length(window_IDs1)
        window1_ID = window_IDs1(i);
        window1_mean_power_spectrum_vector = mean(Zpower(ES_freq_band, window_start_times(window1_ID):window_end_times(window1_ID)), 2);

        for j = 1:length(window_IDs2)
            window2_ID = window_IDs2(j);
            window2_mean_power_spectrum_vector = mean(Zpower(ES_freq_band, window_start_times(window2_ID):window_end_times(window2_ID)), 2);

            % Compute correlation and Fisher z-transform
            r = corr(window1_mean_power_spectrum_vector, window2_mean_power_spectrum_vector, 'type', 'spearman');
            z = 0.5 * log((1 + r) / (1 - r));

            % Store similarity value
            similarity_matrix(i, j) = z;

            % Store mean power spectrum vectors only if needed
            if save_mean_PS_vectors
                window2_mean_PS_vectors(j, :) = window2_mean_power_spectrum_vector;
            end
        end

        % Store mean power spectrum vector for window1 only if needed
        if save_mean_PS_vectors
            window1_mean_PS_vectors(i, :) = window1_mean_power_spectrum_vector;
        end
    end
end
