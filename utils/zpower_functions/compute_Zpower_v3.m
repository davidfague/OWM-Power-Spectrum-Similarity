% v3 adding handling for both across and within trials & use of params.
% v2 revert baselining back to -0.5 - 0 times and across trials
%% note that the first dimension of power is channel and we have changed this to within channel so it is always 1.
function [Zpower] = compute_Zpower_v3(data, params)
    % 1 channel selected; 1 trial selected
    % data should be 1x9001 1xnTimes
    npts = length(data);
    
    % Define parameters for time-frequency analysis
    Fre = params.PSV_freq_band; % Fre = 1:100; % Frequency range from 1 to 100 Hz
    Fsample = params.Fsample; % Fsample = 1000; % Sampling frequency in Hz
    wavenum = params.wavenum; % wavenum = 6; % Number of cycles in the wavelet
    ntrials = size(data,3);

    if (params.baseline_across_trials & (ntrials == 1)) || (~params.baseline_across_trials & (ntrials ~= 1))
        error("The number of trials passed to compute_Zpower_v2 does not match the desired baselining strategy")
    end

    power = zeros(ntrials, 1, length(Fre), npts);
    for ntrial =1:ntrials
        [B, T, F] = BOSC_tf_power(data(:,:,ntrial), Fre, Fsample, wavenum); % Compute time-frequency power
        power(ntrial, 1, :, :) = B; % Store the computed power % note second dim is channel; trials x chan x freq x time
    end
    
    clear B % Clear temporary variable

    % Z-score normalization on baseline
    % Define the baseline period (-0.5 to 0 seconds) % -0.5 to 0 seconds is end of ITI
    T = T - 1; % Adjust time vector to start from 0 (or -1 for T in seconds depending how you look at it)
    wm.powspctrm = power; % Store the power spectrum in a structure
    % bt = T >= -0.5 & T < 0; % Find the indices corresponding to the baseline period (last 500 ms of fixation)
    % bt = T >= -0.75 & T < -0.25; % Find the indices corresponding to the baseline period (middle 250-750 ms of 1 sec fixation)
    bt = T >= params.baseline_T_lims(1) & T < params.baseline_T_lims(2);
    base.powspctrm = power(:,:,:,bt); % Extract the baseline power spectrum % 1 is because there is only 1 channel
    Zpower = zbaseline(wm, base, params.baseline_across_trials); % Compute the z-score normalized power
    Zpower = squeeze(Zpower.powspctrm);

    % %% checking baseline & non-baseline Zpower
    % baseline_subset = Zpower(:,:,:,~bt); %subset by baseline time period (-0.5 to 0)
    % % Step 2: Compute the mean across the 3rd and 4th dimensions
    % % First, compute the mean along the 3rd dimension (if necessary)
    % mean_across_3rd_dim = mean(baseline_subset, 3);
    % % Then compute the mean along the 4th dimension
    % mean_across_3rd_4th_dim = mean(mean_across_3rd_dim, 4);
    % % The result will be a 1x1 matrix, containing the mean of the subsetted data
    % disp(mean_across_3rd_4th_dim);

end