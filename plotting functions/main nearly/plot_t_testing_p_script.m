%% load
patient_id = 201908;
enc_id = 1;
load(sprintf('results/t_testing/p%s_enc%s.mat', num2str(patient_id), num2str(enc_id)))
%% plot a subplot for each image_id (anat labels not visible)
figure;
for image_id = 1:9
    subplot(9,1,image_id)
    % Sort by descending order of final_p_by_channels
    [sorted_p, sort_idx] = sort(p_values_by_chan.EMS(:,image_id), 'descend');
    sorted_anat = anat_by_chan(sort_idx);
    % Plot the sorted data
    bar(sorted_p)
    xticks(1:length(sorted_anat)); % Ensure all ticks are present
    xticklabels(sorted_anat)
    ylabel('p-value') % percent real largest cluster sum t greater than surrogate largest cluster sum t.
    xlabel('Channels')
    title('Sorted p-values by Channel')
end

% Rotate x-axis labels for better visibility
xtickangle(90)

%% one column table

% Initialize variables for the table
image_id = [];
p_value = [];
anat = {};

% Create the table by iterating through the data
for img_id = 1:9
    for chan_idx = 1:length(anat_by_chan)
        image_id = [image_id; img_id];
        p_value = [p_value; p_values_by_chan.EMS(chan_idx, img_id)];
        anat = [anat; anat_by_chan{chan_idx}];
    end
end

% Create the table
p_values_table = table(image_id, p_value, anat, ...
    'VariableNames', {'image_id', 'p_value', 'anat'});

% Sort the table by descending p_value
p_values_table = sortrows(p_values_table, 'p_value', 'descend');

% Display the table
disp(p_values_table);

%% table with image_id columns

% Initialize cell array for the table
sorted_table = cell(length(anat_by_chan), 9);

% Iterate over each image_id
for img_id = 1:9
    % Extract p-values for the current image_id
    p_values = p_values_by_chan.EMS(:, img_id);
    
    % Sort p-values and anatomical labels by descending order of p
    [sorted_p, sort_idx] = sort(p_values, 'descend');
    sorted_anat = anat_by_chan(sort_idx);
    
    % Fill the table column with {anat, p} pairs
    for chan_idx = 1:length(sorted_anat)
        sorted_table{chan_idx, img_id} = {sorted_anat{chan_idx}, int64(sorted_p(chan_idx))};
    end
end

% Convert to a MATLAB table with columns for each image_id
column_names = arrayfun(@(x) sprintf('Image_%d', x), 1:9, 'UniformOutput', false);
sorted_table = cell2table(sorted_table, 'VariableNames', column_names);

% Display the table
disp(sorted_table);

%% 

%% Define patient IDs and enc IDs
patient_ids = [201908, 201903, 201905]; % Add your patient IDs here
enc_id = 1; % Encoding ID

% Initialize combined table
all_p_values_table = table();

% Loop through each patient ID
for pid = patient_ids
    % Load the data for the current patient
    file_name = sprintf('results/t_testing/p%s_enc%s.mat', num2str(pid), num2str(enc_id));
    if isfile(file_name)
        load(file_name);
        
        % Initialize variables for the table
        image_id = [];
        p_value = [];
        anat = {};
        patient_col = []; % Column for patient ID
        chan_id = []; % Column for channel ID

        % Create the table by iterating through the data
        for img_id = 1:9
            for chan_idx = 1:length(anat_by_chan)
                image_id = [image_id; img_id];
                p_value = [p_value; p_values_by_chan.EMS(chan_idx, img_id)];
                anat = [anat; anat_by_chan{chan_idx}];
                patient_col = [patient_col; pid]; % Append patient ID
                chan_id = [chan_id; chan_id_by_chan(chan_idx)]; % Append channel ID
            end
        end

        % Create a table for the current patient
        patient_table = table(image_id, p_value, anat, string(patient_col), chan_id, ...
            'VariableNames', {'image_id', 'p_value', 'anat', 'patient_id', 'chan_id'});
        
        % Sort the table by descending p_value
        patient_table = sortrows(patient_table, 'p_value', 'descend');

        % Append the table to the combined table
        all_p_values_table = [all_p_values_table; patient_table];
    else
        fprintf('File not found for patient %d with enc ID %d\n', pid, enc_id);
    end
end

% Sort the table by descending p_value
all_p_values_table = sortrows(all_p_values_table, 'p_value', 'descend');

% Display the combined table
disp(all_p_values_table);

% Optionally, save the combined table to a file
% writetable(all_p_values_table, 'combined_sorted_p_values_table.csv');
