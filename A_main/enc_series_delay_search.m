use_only_all3_correct = true; % false uses all trials.

addpath('../plotting functions')
output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
mf = matfile(fullfile(output_folder,'201907all3_ES.mat'));
mf_table = matfile(fullfile(output_folder,'201907.mat'));

% rows_without_nans = ~any(any(any(isnan(mf.all3_ES_matrix), 2), 3), 4);
label_table = mf_table.label_table;
rows_without_nan = get_rows_without_nan(label_table);
if use_only_all3_correct
    rows_with_all3_correct = sum(label_table.encoding_correctness, 2)==3;
    rows_to_use = rows_without_nan & rows_with_all3_correct;
else
    rows_to_use = rows_without_nan;
end

label_table = label_table(rows_to_use,:);

indices_to_use = find(rows_to_use);
max_id = max(indices_to_use);
min_id = min(indices_to_use);

colors = {'r', 'g', 'y'}; % enc1, enc2, enc3 colors
enc_labels = {'Encoding 1 mean', 'Encoding 1 std','Encoding 2 mean', 'Encoding 2 std', 'Encoding 3 mean', 'Encoding 3 std'};

fig = figure();
hold on
for enc_id = 1:3

    ES_matrix = mf.all3_ES_matrix(min_id:max_id, :,:,enc_id);
    ES_matrix = ES_matrix(indices_to_use - (min_id-1),:,:,:);
    
    % flatten across encoding
    % size(ES_matrix)
    % ES_matrix = mean(ES_matrix, 2);
    plot_avg_ES_over_time(ES_matrix, fig, colors{enc_id}, enc_labels{enc_id}, use_only_all3_correct);
end

legend(enc_labels, 'Location', 'best'); % Add legend

if use_only_all3_correct
    savefig(fig, '../results/ES_avg_trace_all3correct_trials.fig')
else
    savefig(fig, '../results/ES_avg_trace_all_trials.fig')
end
hold off;

function plot_avg_ES_over_time(ES_matrix, fig, color, label, use_only_all3_correct)
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
    figure(fig);
    hold on

    % Plot average line
    plot(ntime, avg_ES_over_time, 'Color', color, 'LineWidth', 2, 'DisplayName', label);

    % Plot shaded area for standard deviation
    fill([ntime, fliplr(ntime)], [upper_bound', fliplr(lower_bound')], ...
         color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Customize plot
    xlabel('Time Window'); % change to time
    ylabel('Encoding Similarity');
    if use_only_all3_correct
        title('Encoding similarity during all 3 correct trials');
    else
        title('Encoding similarity during all trials');
    end
    grid on;
    add_event_lines()
    options.DashedLines = true;
    options.LineColor = color; % Red color
    plot_horizontal_means(avg_ES_over_time, options)
    % hold off;
end
