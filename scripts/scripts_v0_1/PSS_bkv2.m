%% calc z score power all patients

clear all
close all

addpath '/Users/Bornali/Library/CloudStorage/OneDrive-UniversityofMissouri/Documents - Kundu Lab - Ogrp/data/RSA_analysis/Code/Power Spectrum Similarity/subfunctions'


str1='/Users/Bornali/Library/CloudStorage/OneDrive-UniversityofMissouri/Documents - Kundu Lab - Ogrp/data/OWM Utah Data/';
ptfiles={  'CS201901';
    'CS201902';
    'CS201903';
    'CS201905';
    'CS201906';
    'CS201907';
    'CS201908';
    'CS201910';
    'CS201915';
    'UIC202117'};

%%

for pt=1%:length(ptfiles)
    clearvars -except pt ptfiles usebi str1
    dirname=[str1 ptfiles{pt}];
    cd(dirname)
    %mkdir('zpower')

    load('D_OWM_t_bipolar.mat')
    % load('OWM_trialinfo.mat', 'Cond_performance')
    load('D_OWM_t_bipolar.mat', 'labelsanatbkedit')
    load('gammachans2sd_alltrials.mat', 'sigchans2')

    num_channels = size(D_OWM_t,1);
    num_trials = size(D_OWM_t,3);
    num_timepoints = size(D_OWM_t,2); % 9001 time points

    Zpower_all=nan( num_trials, 40, num_timepoints);

    for channel_idx = 1:num_channels
        Zpower_all=nan( num_trials, 40, num_timepoints);
        for trial_idx = 1:num_trials
            data_subset=D_OWM_t(channel_idx, :,trial_idx);
            if all(isnan(data_subset))
                Zpower=nan(40,num_timepoints);
            else
                Zpower = compute_Zpower(data_subset);

            end
            Zpower_all( trial_idx,:,:)=Zpower;

            clear Zpower data_subset
        end
        zfile_chan=join([labelsanatbkedit.ROI{channel_idx}, 'zpow.mat'], '_')
        save(join([dirname, '/zpower/', zfile_chan], ''), 'Zpower_all')
    end

end

%% Similarity processing calc ES for sig chans 2
clearvars -except  ptfiles str1

ES_prefix = 'wholeTrial';
time = 1:9001;
% 100ms time window with step of 10ms
wt = length(time);
num_windows = floor((wt - 100) / 10) + 1; % Correct calculation for number of windows
% Precompute window start and end times
window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
clear num_windows wt
window_end_times = time(window_start_times + 99);            % End times

% access window start and end: (:, window_id)
% windows for 2nd iteration (whole trial windows)
window_IDs2 = true(size(time));
window_IDs2 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs2, window_start_times, window_end_times);
% Define the time windows for each encoding period
stim1_start = 1000; stim1_end = 1500; % stim 1 time window
stim2_start = 1500; stim2_end = 2000; % stim 2 time window
stim3_start = 2000; stim3_end = 2500; % stim 3 time window

ES_freq_band = 1:40;
save_mean_PS_vectors = false;
mkdir('similarity')
for pt=1%:length(ptfiles)
    clearvars -except pt ptfiles str1
    dirname=[str1 ptfiles{pt}];
    cd(dirname)
    %mkdir('zpower')

    % load('D_OWM_t_bipolar.mat')
    % load('OWM_trialinfo.mat', 'Cond_performance')
    load('D_OWM_t_bipolar.mat', 'labelsanatbkedit')
    load('gammachans2sd_alltrials.mat', 'sigchans2')

    num_channels = size(labelsanatbkedit,1);

    for channel_idx = 1:num_channels
        if any(channel_idx==sigchans2)

            zfile_chan=join([labelsanatbkedit.ROI{channel_idx}, 'zpow.mat'], '_')
            load(join([dirname, '/zpower/', zfile_chan], ''), 'Zpower_all')

            num_trials = size(Zpower_all,1);
            num_timepoints = size(Zpower_all,3); % 9001 time points
       
            for enc_ID=1:3
                % get ntimes logical array corresponding to this
                % encoding timeframe
                window_IDs1 = get_encoding_times_from_enc_ID(enc_ID, time, stim1_start, stim1_end, stim2_start, stim2_end, stim3_start, stim3_end);
                % get window IDs from logical array
                window_IDs1 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs1, window_start_times, window_end_times);
                % now have 2 logical arrays of size nTimes:
                % window_IDs1, window_IDs2
                %use logical arrays to subset window_IDs.
                sim_all=nan(num_trials,3, 40,891);
                for tr=1:num_trials
                    [similarity_matrix, window1_mean_PS_vectors, window2_mean_PS_vectors] = compute_ES(squeeze(Zpower_all(tr,:,:)), window_start_times, window_end_times, window_IDs1, window_IDs2, ES_freq_band, save_mean_PS_vectors);
                    sim_all(tr,enc_ID,:,:)=similarity_matrix;
                    clear similarity_matrix
                end

            end
           
            ES_file_chan=join([labelsanatbkedit.ROI{channel_idx}, 'sim.mat'], '_')

            save(join([dirname, '/similarity/', ES_file_chan], ''), "sim_all",  '-v7.3')

        end
    end

end



