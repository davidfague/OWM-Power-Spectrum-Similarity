% look at what channels can discern multiple images
% will need to add session id to patient_channel_str for missouri data
table_to_load = 'summary_p_tables_midbase_1-40hz_utah_clip_infs.mat';
load(table_to_load)

patient_channel_str = strings(length(all_p_table.patient_id), 1);

for i =1:length(all_p_table.patient_id)
    if all_p_table.chan_id(i) < 10
        patient_channel_str(i) = sprintf("p%s_ch0%s", ...
            num2str(all_p_table.patient_id(i)), num2str(all_p_table.chan_id(i)));
    else
        patient_channel_str(i) = sprintf("p%s_ch%s", ...
            num2str(all_p_table.patient_id(i)), num2str(all_p_table.chan_id(i)));
    end
end

all_p_table.patient_channel_str = patient_channel_str;
clear patient_channel_str

unique(all_p_table.patient_channel_str);

significant_table = all_p_table(all_p_table.p>0.975,:);

unique_patient_channel_strs = unique(significant_table.patient_channel_str);

anats = strings(length(unique_patient_channel_strs), 1);
numbers = nan(length(unique_patient_channel_strs), 1);
% sig_images = cells(nan(length(unique_patient_channel_strs), 1));

for pat_chan_idx = 1:length(unique_patient_channel_strs)
    current_pat_chan = unique_patient_channel_strs(pat_chan_idx);
    significant_table_current_chan = significant_table(strcmp(current_pat_chan, significant_table.patient_channel_str),:);
    anat = unique(significant_table_current_chan.anat);
    anats(pat_chan_idx) = anat;
    number_significant_images_current_chan = size(significant_table_current_chan, 1);
    numbers(pat_chan_idx) = number_significant_images_current_chan;
    % sig_images = unique(significant_table_current_chan.image_id)';
end

num_sig_images_by_chan_table = table(unique_patient_channel_strs, anats, numbers, ...
    'VariableNames',{'patient_channel_str', 'anat', 'num_sig_images'});

save(table_to_load, "num_sig_images_by_chan_table", "-append")

%% filter table to keep channels that differentiate > 1 image or work with neighboring channels to differentiate > 1 image
% 1. Parse the patient ID and channel number from patient_channel_str.
%    Assume strings are in the format 'pYYYYMM_chNN'
tokens = regexp(num_sig_images_by_chan_table.patient_channel_str, '^(p\d+)_ch(\d+)$', 'tokens');
patientIDs = cellfun(@(x) x{1}{1}, tokens, 'UniformOutput', false);
channelNums = cellfun(@(x) str2double(x{1}{2}), tokens);

% 2. Add new columns for easier processing
num_sig_images_by_chan_table.patientID = patientIDs;
num_sig_images_by_chan_table.channelNum = channelNums;

% 3. Initialize a logical index: keep rows with >2 num_sig_images
keepFlag = num_sig_images_by_chan_table.num_sig_images > 2;

% 4. For rows that might be kept based on being "close to neighbors", group by patientID and anat.
[G, ~] = findgroups(num_sig_images_by_chan_table.patientID, num_sig_images_by_chan_table.anat);

% Loop over each group and flag rows that are within 5 channels of a neighbor.
uniqueGroups = unique(G);
for g = uniqueGroups'
    groupIdx = find(G == g);
    if numel(groupIdx) > 1
        % Sort the group by channel number.
        [sortedCh, sortOrder] = sort(num_sig_images_by_chan_table.channelNum(groupIdx));
        
        % Compute differences between successive channels.
        diffCh = diff(sortedCh);
        
        % Create a logical array to flag adjacent channels within 5.
        % If the difference between adjacent channels is <= 5, mark both as "close".
        closeFlag = diffCh <= 5;
        groupKeep = false(size(sortedCh));
        groupKeep(1:end-1) = groupKeep(1:end-1) | closeFlag;
        groupKeep(2:end)   = groupKeep(2:end)   | closeFlag;
        
        % Map the flag back to the original table indices.
        % sortOrder gives the order within groupIdx.
        keepFlag(groupIdx(sortOrder(groupKeep))) = true;
    end
end

% 5. Finally, filter the table.
many_sig_images_by_chan_table = num_sig_images_by_chan_table(keepFlag, :);

% get rid of rows that have anat = "Unknown", "Right Cerebral White Matter"
% "Resection cavity" "Left Cerebral White Matter"
stringsToRemove = {'Unknown', 'Right Cerebral White Matter', 'Resection cavity', 'Left Cerebral White Matter'};
rowsToRemove = ismember(many_sig_images_by_chan_table.anat, stringsToRemove);
many_sig_images_by_chan_table(rowsToRemove, :) = [];

% fix patient ID
many_sig_images_by_chan_table.patientID = double(erase(string(many_sig_images_by_chan_table.patientID), "p"));

% back track and find the significant images for the outstanding channels
rows_to_keep = ismember(significant_table.patient_channel_str, many_sig_images_by_chan_table.patient_channel_str);
outstanding_obs_table = significant_table(rows_to_keep,:);

save(table_to_load, "many_sig_images_by_chan_table", "outstanding_obs_table", "-append")