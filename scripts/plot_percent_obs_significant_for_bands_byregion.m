% this version will adjust the npatient_per_region filtering based on the
% data subset

%% parameters
analysis_params = struct();
analysis_params.min_npatient_per_region = 2;
% Set to true to filter regions by n_unique_patients > min_npatient_per_region;
% set to false to use the predefined regions from regions_to_use.mat.
analysis_params.use_nunique_filter = true;

p_threshhold = 0.975;
threshold_str = sprintf('p > %.3f', p_threshhold);

% Define frequency bands to analyze
% es_freq_bands = {1:40, 8:12, 12:30, 30:70};  % each cell is a vector of frequencies
es_freq_bands = {1:40, 1:8, 8:20, 20:40};

% Data set selection
k12wm = false;
if k12wm
    data_set_str = 'k12wm';
else
    data_set_str = 'utah';
end

custom_params = struct();
custom_params.hellbender = true;
params = get_parameters(custom_params);

%% Load regions_to_use if not filtering by unique patients
if ~analysis_params.use_nunique_filter
    if params.hellbender
        load("../results/regions_to_use.mat", "regions_to_use")
    else
        load("..\results\regions_to_use.mat", "regions_to_use")
    end
end

%% Process each frequency band and compute summary statistics
numBands = numel(es_freq_bands);
summaryTables = cell(numBands, 1);
regionLists = cell(numBands, 1); % store filtered region names for each band

for fb = 1:numBands
    % Define current frequency band and string for file naming
    current_band = es_freq_bands{fb};
    freq_min = min(current_band);
    freq_max = max(current_band);
    freq_str = sprintf("%d-%dhz", freq_min, freq_max);
    
    % Construct file name and load the data table
    table_to_load = sprintf('summary_p_tables_midbase_%s_%s_clip_infs.mat', freq_str, data_set_str);
    loaded_table = load(table_to_load);
    all_patients_p_table = loaded_table.all_p_table;
    
    %% Preprocessing: remove missing and standardize labels
    all_patients_p_table = all_patients_p_table(~ismissing(all_patients_p_table.anat), :);
    all_patients_p_table.anat = lower(all_patients_p_table.anat);
    all_patients_p_table.anat_merged = removeDirectionPrefix(all_patients_p_table.anat);
    
    %% Compute summary statistics per anatomical region
    regions_in_table = unique(all_patients_p_table.anat_merged);
    nRegions_table = numel(regions_in_table);
    
    totalOccurrences    = zeros(nRegions_table, 1);
    fracAboveThreshold  = zeros(nRegions_table, 1);
    totalAboveThreshold = zeros(nRegions_table, 1);
    nUniquePatients     = zeros(nRegions_table, 1);
    
    for i = 1:nRegions_table
        idx = strcmp(all_patients_p_table.anat_merged, regions_in_table{i});
        totalOccurrences(i)    = sum(idx);
        fracAboveThreshold(i)  = sum(all_patients_p_table.p(idx) > p_threshhold) / totalOccurrences(i);
        totalAboveThreshold(i) = sum(all_patients_p_table.p(idx) > p_threshhold);
        nUniquePatients(i)     = numel(unique(all_patients_p_table.patient_id(idx)));
    end
    
    % Create a summary table for the current frequency band
    summaryTable = table(regions_in_table, totalOccurrences, totalAboveThreshold, fracAboveThreshold, nUniquePatients, ...
        'VariableNames', {'anat_merged', 'total_occurrences', 'total_above_threshold', 'frac_above_threshold', 'n_unique_patients'});
    
    % Apply filtering based on user choice
    if analysis_params.use_nunique_filter
        % Filter regions by n_unique_patients > threshold
        summaryTable = summaryTable(summaryTable.n_unique_patients > analysis_params.min_npatient_per_region, :);
    else
        % Use regions_to_use filtering
        summaryTable = summaryTable(ismember(summaryTable.anat_merged, regions_to_use), :);
        % Add placeholders for any missing regions
        for region_idx = 1:length(regions_to_use)
            region = regions_to_use{region_idx};
            if ~ismember(region, summaryTable.anat_merged)
                summaryTable(end+1,:) = table({region}, 0, nan, nan, 0, ...
                    'VariableNames', {'anat_merged', 'total_occurrences', 'total_above_threshold', 'frac_above_threshold', 'n_unique_patients'});
            end
        end
        % Sort the table to match the order in regions_to_use
        [~, sortIdx] = ismember(regions_to_use, summaryTable.anat_merged);
        summaryTable = summaryTable(sortIdx, :);
    end
    
    % Store the summary table and region list for this frequency band
    summaryTables{fb} = summaryTable;
    regionLists{fb} = summaryTable.anat_merged;
end

%% Determine regions to plot across frequency bands
if analysis_params.use_nunique_filter
    % Use intersection of regions across all frequency bands
    common_regions = regionLists{1};
    for fb = 2:numBands
        common_regions = intersect(common_regions, regionLists{fb});
    end
    regions_to_plot = common_regions;
else
    regions_to_plot = regions_to_use;
end

%% Combine summary statistics across frequency bands for plotting
numRegions = numel(regions_to_plot);
frac_matrix = nan(numRegions, numBands);  % rows: regions; columns: frequency bands

for fb = 1:numBands
    st = summaryTables{fb};
    % For each region in regions_to_plot, fill in the fraction value if available
    for r = 1:numRegions
        idx = strcmp(st.anat_merged, regions_to_plot{r});
        if any(idx)
            frac_matrix(r, fb) = st.frac_above_threshold(idx);
        end
    end
end

%% Plot grouped bar chart (vertical)
figure('WindowState', 'maximized');
b = bar(frac_matrix, 'grouped');  % groups: regions; bars: frequency bands
set(gca, 'XTick', 1:numRegions, 'XTickLabel', regions_to_plot, 'FontSize', 12);
xlabel('Anatomical Region');
ylabel(sprintf('Fraction of Significant Observations %s', threshold_str));
title_str = sprintf('Temporally Generalized Significant WI-BI\n%s dataset', data_set_str);
if analysis_params.use_nunique_filter
    title_str = sprintf('%s\n(Filtered: n_{unique patients} > %d)', title_str, analysis_params.min_npatient_per_region);
end
title(title_str);

% Create legend entries using frequency band labels
legendLabels = cellfun(@(band) sprintf('%d-%dhz', min(band), max(band)), es_freq_bands, 'UniformOutput', false);
legend(legendLabels, 'Location', 'bestoutside');
