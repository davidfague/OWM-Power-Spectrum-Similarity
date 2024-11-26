%% make sure all 3 correct
% plot trial ES for all 3 enc
% for trials with all 3 correct
addpath('../plotting functions')
output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
mf = matfile(fullfile(output_folder,'201907all3_ES.mat'));
mf_table = matfile(fullfile(output_folder,'201907.mat'));

% rows_without_nans = ~any(any(any(isnan(mf.all3_ES_matrix), 2), 3), 4);
label_table = mf_table.label_table;
rows_without_nan = get_rows_without_nan(label_table);

label_table = label_table(rows_without_nan,:);

indices_without_nan = find(rows_without_nan);
max_id = max(indices_without_nan);
min_id = min(indices_without_nan);

fig = figure();
hold on
for enc_id = 1:3

    ES_matrix = mf.all3_ES_matrix(min_id:max_id, :,:,enc_id);
    ES_matrix = ES_matrix(indices_without_nan - (min_id-1),:,:,:);
    
    % flatten across encoding
    % size(ES_matrix)
    % ES_matrix = mean(ES_matrix, 2);
    plot_avg_ES_over_time(ES_matrix, fig)
end

function plot_avg_ES_over_time(ES_matrix, figure)
    avg_ES_over_time = squeeze(mean(ES_matrix, [1 2]));
    % avg_ES_over_time(isinf(avg_ES_over_time)) = 0;
    ntime = 1:length(avg_ES_over_time); % can update to do the mean of window start and end time for this windowID.
    
    std_ES_over_time = squeeze(std(ES_matrix, 0, [1 2]));
    clear ES_matrix
    std_ES_over_time(isnan(std_ES_over_time)) = 0;
    
    % figure()
    % figure.plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
    y_limits = ylim;
    avg_ES_over_time(isinf(avg_ES_over_time)) = y_limits(2);
    
    upper_bound = avg_ES_over_time + std_ES_over_time;
    lower_bound = avg_ES_over_time - std_ES_over_time;
    
    figure.plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
    hold on;
    
    % Plot the shaded area using the fill function
    figure.fill([ntime, fliplr(ntime)], [upper_bound', fliplr(lower_bound')], ...
         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Customize plot
    xlabel('Time Window'); % change to time
    ylabel('Encoding Similarity');
    title('Encoding similarity during all trials');
    grid on;
    add_event_lines()
    plot_horizontal_means(avg_ES_over_time)
    % hold off;
end
