function plot_and_save_RSA(PS_file, output_folder)
    % PLOT_AND_SAVE_RSA - This function filters, processes, and saves results for Power Spectrum Similarity RSA analysis
    
    % Step 1: Extract the label table and remove NaN rows
    label_table = PS_file.label_table;
    rows_with_nan = any(isnan(label_table.EMS_means), 2);
    label_table = label_table(~rows_with_nan, :);
    
    % Step 2: Subset the table based on encoding performance
    label_table = subset_table_by_enc_performance(label_table);
    
    % Step 3: Extract PSVs and plot whole-trial PSVs for each anatomy and region
    % Get all windowed mean PS vectors and filter by regions/channels/trials as necessary
    rows_without_nan = find(~rows_with_nan);
    all_PS_vectors = PS_file.all_windowed_mean_PS_vectors(:, :, :);
    all_PS_vectors = all_PS_vectors(:, :, rows_without_nan);
    
    % Step 4: Plot PSVs by anatomical regions and save plots
    anatomical_labels = unique(string(label_table.anatomical_label(:)));
    for i = 1:length(anatomical_labels)
        anat_label = anatomical_labels{i};
        
        % Filter by anatomical region
        region_indices = strcmp(label_table.anatomical_label, anat_label);
        region_PSVs = all_PS_vectors(:, :, region_indices);
        
        % Plot mean PSVs for this region
        figure;
        plot(mean(region_PSVs, 3, 'omitnan'));
        title(['Whole-Trial PSVs for ', anat_label]);
        xlabel('Time Points');
        ylabel('Power Spectrum');
        
        % Save the plot
        plot_filename = fullfile(output_folder, ['PSV_', anat_label, '.png']);
        fprintf('Finish %s', plot_filename)
        saveas(gcf, plot_filename);
        close(gcf);
    end
    
    % Step 5: Calculate EMS and EFS for the cleaned data
    [EMS, EFS] = compute_EMS_EFS(PS_file, all_PS_vectors);
    
    % Step 6: Save results to the output folder
    results_file = fullfile(output_folder, ['RSA_Results_', num2str(PS_file.patient_ID), '.mat']);
    save(results_file, 'EMS', 'EFS', 'label_table', 'anatomical_labels');
end