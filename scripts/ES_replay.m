%% For computing Phase-Amplitude Coupling of Encoding similarity, that is amplitude is encoding similarity.
% Phase-Similarity Coupling
% assuming in AnalyzingNeuralTimeSeries 
% and AnalyzingNeuralTimeSeries and OWM-Power-Spectrum-Similarity are in
% the same folder
% addpath("../OWM-Power-Spectrum-Similarity/scripts/")

% actually should be in "/OWM-Power-Spectrum-Similarity/scripts/" for
% getting params and data.

%% parameters

script_specs = struct();
script_specs.clip_inf_similarities = true;
script_specs.skip_existing = false;

custom_params = get_custom_params(script_specs); % get other defaults like hellbender, output_folder_name; (intending to .gitignore 'get_custom_params.m')

params = get_parameters(custom_params);

clear script_specs custom_params

%%

params.patient_id = 201910;
params.chan_id = 23;
params.image_id = 3;
params.enc_window_ids = params.enc1_win_IDs; % change with enc_id
params.enc_id = 1;
params.session_id = 1;

params.type = 'EMS'; % for plotting WI_vs_BI only
params.mean_out_time_dimensions = true;
params.average_diff = true;
params.same_n = true; % only affects plot_similarity_means_heatmaps

params.only_all3_correct = true; % filters test trials (with item)

band_idx = 2;%:length(params.bands)
params.band_to_process = params.freq_band_map(params.bands{band_idx});
params.freq_min = min(params.band_to_process.range);
params.freq_max = max(params.band_to_process.range);

%% params that will need updated if their input is changed below
% might need to put the above params in script_specs and then put this
% logic toward end of get_parameters and use 'if the required fields are present'.

params.patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, params.patient_id);
params.patient_preprocessed_data_path = params.patient_preprocessed_data_paths{params.session_id};

params.WI_BI_folder_to_save_in = fullfile(sprintf("results/WI vs BI/p%s chan%s image%s enc%s sess%d", ...
            num2str(params.patient_id), num2str(params.chan_id), ...
            num2str(params.image_id), num2str(params.enc_id), params.session_id));

params.anat_labels = get_anat_labels(params.patient_preprocessed_data_path, params);

% params.image_labels = load(fullfile(params.patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted

params.anat = string(params.anat_labels.labelsanatbkedit.anatmacro1(params.chan_id));

%% load signal

[WI, BI] = load_WI_BI_for_channel(params);

[WI, BI] = clip_infs_of_z_similarities(WI, BI);

data = mean(WI,1);
% 
% trial_type = 'single_trial'; % can be either single_trial or mean_trial, can loop through both.
% trial_idx = 10; % trial pair for single trial
% 
% if strcmp(trial_type, 'single_trial')
%     corr_vector = data(:,:,trial_idx); % pick one trial pair
%     corr_vector = squeeze(mean(corr_vector, 1)); % mean across encoding times
% elseif strcmp(trial_type, 'mean_trial')
%     corr_vector = squeeze(mean(data, [1 3])); % mean across trial pairs and encoding times
% end

%% structure as signal as assumed

types = {'mean_across_trials', 'single_trial', 'all_trials'};
analysis_type = types{1};
single_trial_index = 1;

EEG = struct();
EEG.srate = 1000/10; % 10 ms windows;
EEG.Fs = EEG.srate;
EEG.data(1,:,:) = data; % channels, times, trials
EEG.pnts = size(EEG.data,2);
EEG.npoints = EEG.pnts;
EEG.nbchan = size(EEG.data,1);
EEG.trials = size(EEG.data,3);

EEG.times = 2.6:1/EEG.srate:6.4; % maintenance period
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
if strcmp(analysis_type, 'all_trials') % adjust based on trials
    n_data        = EEG.pnts*EEG.trials; % multiple trials
else
    n_data        = EEG.pnts*1;
end
n_convolution = n_wavelet+n_data-1; 

% FFT of data
% fft_EEG = fft(reshape(EEG.data(1,:,:),1,EEG.pnts*EEG.trials),n_convolution);
if strcmp(analysis_type, 'all_trials')
    EEG.signal_to_fft = EEG.data(1,:,:);
    EEG.fft = fft(reshape(signal_to_fft,1,EEG.pnts*EEG.trials),n_convolution);
elseif strcmp(analysis_type, 'mean_across_trials')
    EEG.signal_to_fft = squeeze(mean(EEG.data(1,:,:), 3));
    EEG.fft = fft(signal_to_fft, n_convolution);
elseif strcmp(analysis_type, 'single_trial')
    EEG.signal_to_fft = squeeze(EEG.data(1,:,single_trial_index));
    EEG.fft = fft(signal_to_fft, n_convolution);
else
    error('analysis_type %s Not Implemented', analysis_type)
end

in_dB = true;
title_prefix = sprintf("%s %s p%s chan%s image%s enc%s\n %s", ...
        params.type, analysis_type, num2str(params.patient_id), ... 
        num2str(params.chan_id), num2str(params.image_id), ...
        num2str(params.enc_id), params.anat);
windowed = false;

fig = plot_fft(EEG, in_dB, windowed, title_prefix);

% vars_to_keep = whos().name; % if wanting to keep all variables from up till now.
%% Figure 1, plotting the EMS signal, EMS power of a frequency
% Two options:
% (1) compute PAC between the phase of some frequency and the amplitude of
% the EMS itself
% (2) compute PAC between the phase of some frequency and the amplitude of
% the power of some frequency of the EMS.

%% option 2
freq4phase = 5; % in Hz
freq4power = 5; 

% wavelet for phase and its FFT
wavelet4phase = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
fft_wavelet4phase = fft(wavelet4phase,n_convolution);

% wavelet for power and its FFT
wavelet4power = exp(2*1i*pi*freq4power.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4power))^2));
fft_wavelet4power = fft(wavelet4power,n_convolution);

% get phase values
convolution_result_fft = ifft(fft_wavelet4phase.*fft_EEG,n_convolution);
phase = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));
if all(phase == 0) || any(isnan(phase))
    error("")
end
% get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
convolution_result_fft = ifft(fft_wavelet4power.*fft_EEG,n_convolution);
pwr = abs(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size)).^2;

% plot power and phase
figure
subplot(131)
plot(EEG.times(1:end-1),phase)
hold on
plot(EEG.times(1:end-1),(pwr-mean(pwr))/std(pwr),'r')
legend({[num2str(freq4phase),' Hz phase'];[num2str(freq4power),' Hz power']})
set(gca,'xlim',[min(EEG.times) max(EEG.times)])
axis square

% plot power as a function of phase in polar space
subplot(132)
polarplot(phase,pwr,'.')

% plot histogram of power over phase
n_hist_bins = 30;

phase_edges=linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);

for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(pwr(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

subplot(133)
bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)
axis square

%% option 1
freq4phase = 15; % in Hz
freq4power = 'EMS'; 
pwr = signal_to_fft;

% wavelet for phase and its FFT
wavelet4phase = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
fft_wavelet4phase = fft(wavelet4phase,n_convolution);

% % wavelet for power and its FFT
% wavelet4power = exp(2*1i*pi*freq4power.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4power))^2));
% fft_wavelet4power = fft(wavelet4power,n_convolution);

% get phase values
convolution_result_fft = ifft(fft_wavelet4phase.*fft_EEG,n_convolution);
phase = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));
if all(phase == 0) || any(isnan(phase))
    error("")
end
% get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
% convolution_result_fft = ifft(fft_wavelet4power.*fft_EEG,n_convolution);
% pwr = abs(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size)).^2;

% plot power and phase
figure
subplot(131)
plot(EEG.times(1:end-1),phase)
hold on
plot(EEG.times(1:end-1),(pwr-mean(pwr))/std(pwr),'r')
legend({[num2str(freq4phase),' Hz phase'];freq4power})
set(gca,'xlim',[min(EEG.times) max(EEG.times)])
axis square

% plot power as a function of phase in polar space
subplot(132)
polarplot(phase,pwr,'.')

% plot histogram of power over phase
n_hist_bins = 30;

phase_edges=linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);

for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(pwr(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

subplot(133)
bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Mean ' num2str(freq4power) '' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)
axis square

%% measure PAC and PACz
%% Figure 30.4 (intentionally created a biased signal by removing the 1/4 phase)

phase_bias = phase;
power_bias = pwr;

phase_bias(phase<-pi/2) = [];
power_bias(phase<-pi/2) = [];

% figure
% % plot power as a function of phase in polar space
% subplot(131)
% polar(phase,pwr,'.')
% title([ 'PAC = ' num2str(round(abs(mean(pwr.*exp(1i*phase))))) ])
% 
% subplot(132)
% polar(phase,pwr*10,'.')
% title([ 'PAC = ' num2str(round(abs(mean(pwr*10.*exp(1i*phase))))) ])
% 
% subplot(133)
% polar(phase_bias,power_bias,'.')
% title([ 'PAC = ' num2str(round(abs(mean(power_bias.*exp(1i*phase_bias))))) ])
%% Figure 30.5 (showing the PACz for the data and the intentionally biased data)

% observed cross-frequency-coupling (note the similarity to Euler's formula)
obsPAC = abs(mean(pwr.*exp(1i*phase)));
obsPAC_bias = abs(mean(power_bias.*exp(1i*phase_bias)));

num_iter = 1000;

permutedPAC = zeros(2,num_iter);

for i=1:num_iter
    
    % select random time point
    random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
    random_timepoint_bias = randsample(round(length(power_bias)*.8),1)+round(length(power_bias)*.1);
    
    % shuffle power
    timeshiftedpwr      = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
    timeshiftedpwr_bias = [ power_bias(random_timepoint_bias:end) power_bias(1:random_timepoint_bias-1) ];
    
    % compute PAC
    permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
    permutedPAC(2,i) = abs(mean(timeshiftedpwr_bias.*exp(1i*phase_bias)));
end

% compute PACz
pacz(1) = (obsPAC-mean(permutedPAC(1,:)))/std(permutedPAC(1,:));
pacz(2) = (obsPAC_bias-mean(permutedPAC(2,:)))/std(permutedPAC(2,:));

figure
subplot(221)
hist(permutedPAC(1,:),50);
hold on
plot([obsPAC obsPAC],get(gca,'ylim')/2,'m','linew',3)
legend({'Permuted values';'Observed value'})
xlabel('Modulation strength'), ylabel('Number of observations')
title([ 'PAC_z = ' num2str(pacz(1)) ])

subplot(222)
hist(permutedPAC(2,:),50)
hold on
plot([obsPAC_bias obsPAC_bias],get(gca,'ylim')/2,'m','linew',3)
legend({'Permuted values';'Observed value'})
xlabel('Modulation strength'), ylabel('Number of observations')
title([ 'PAC_z = ' num2str(pacz(2)) 'intentionally biased by removing 1/4 phase.'])

% plot histogram of power over phase
n_hist_bins = 30;
phase_edges=linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);
for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(pwr(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

subplot(223)
h=bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(h,'linestyle','none'); % turn off black lines around histogram bars
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)


% plot histogram of power over phase
n_hist_bins = 30;
phase_edges=linspace(min(phase_bias),max(phase_bias),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);
for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(power_bias(phase_bias>phase_edges(i) & phase_bias<phase_edges(i+1)));
end

subplot(224)
h=bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(h,'linestyle','none'); % turn off black lines around histogram bars
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)

%% Figure 30.6
num_iter = 1000;
permutedPAC = zeros(3,num_iter);
eeg = signal_to_fft;
pwr = signal_to_fft;
obsPAC = abs(mean(pwr.*exp(1i*phase)));
for i=1:num_iter
    
    % Permutation method 1: select random time point
    random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
    timeshiftedpwr   = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
    permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
    
    % Permutation method 2: totally randomize power time series
    permutedPAC(2,i) = abs(mean(pwr(randperm(length(pwr))).*exp(1i*phase)));
    
    % Permutation method 3: FFT-based power time series randomization
    f = fft(pwr); % compute FFT
    A = abs(f);   % extract amplitudes
    zphs=cos(angle(f))+1i*sin(angle(f)); % extract phases
    powernew=real(ifft(A.*zphs(randperm(length(zphs))))); % recombine using randomized phases (note: use original phases to prove that this method reconstructs the original signal)
    powernew=powernew-min(powernew);
    
    permutedPAC(3,i) = abs(mean(powernew.*exp(1i*phase)));
end

% compute PACz and plot
figure
for i=1:3
    subplot(2,3,i)
    
    % plot example power time series
    switch i
        case 1
            plot(EEG.times(1:end-1),timeshiftedpwr)
            title('H_0: Time-shifted')
        case 2
            plot(EEG.times(1:end-1),pwr(randperm(length(pwr))))
            title('H_0: randomized')
        case 3
            plot(EEG.times(1:end-1),powernew)
            title('H_0: FFT-derived randomization')
    end
    set(gca,'xlim',[0 EEG.times(end)],'ylim',[min(pwr) max(pwr)])
    ylabel("EMS")
    
    % plot null-hypothesis distribution
    subplot(2,3,i+3)
    pacz = (obsPAC-mean(permutedPAC(i,:)))/std(permutedPAC(i,:));
    [y,x]=hist(permutedPAC(i,:),50);
    h=bar(x,y,'histc');
    set(h,'linestyle','none');
    hold on
    plot([obsPAC obsPAC],get(gca,'ylim')/2,'m','linew',3)
    legend({'Permuted values';'Observed value'})
    xlabel('Modulation strength'), ylabel('Number of observations')
    title([ 'PAC_z = ' num2str(pacz) ])
end

sgtitle(['Methods for generating null distribution. Null Hypothesis: There is no temporal relationship between phase at ' num2str(freq4phase) ' Hz and EMS'])


%% measure PACz across multiple frequencies

EEG = struct();
EEG.srate = 1000/10; % 10 ms windows;
EEG.data(1,:,:) = data; % channels, times, trials
EEG.pnts = size(EEG.data,2);
EEG.nbchan = size(EEG.data,1);
EEG.trials = size(EEG.data,3);

EEG.times = 2.6:1/EEG.srate:6.4; % maintenance period
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*1;
n_convolution = n_wavelet+n_data-1; 
num_frequencies = 1000;
frequencies = get_frequencies(num_frequencies, [1e-15 30], 'log');
PACz = nan(num_frequencies);

apply_window = false;

freq4power = 'EMS'; 

signal_to_fft = squeeze(mean(EEG.data(1,:,:), 3));

if apply_window
    N = length(signal_to_fft);              % Length of the signal
    w = hann(N)';               % Generate a Hann window (note the transpose to match dimensions
    % Apply the window to the signal
    signal_to_fft = signal_to_fft .* w;
end

fft_EEG = fft(signal_to_fft, n_convolution);

num_iter = 1000;
permutedPAC = zeros(3,num_iter);
pacz = nan(3, num_frequencies);

eeg = signal_to_fft;
pwr = signal_to_fft;

for freq_idx = 1:num_frequencies
    freq4phase = frequencies(freq_idx); % in Hz
    
    % calculate phase
    % wavelet for phase and its FFT
    wavelet4phase = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
    fft_wavelet4phase = fft(wavelet4phase,n_convolution);

    %convolve and ifft phase's wavelet with signal's fft
    convolution_result_fft = ifft(fft_wavelet4phase.*fft_EEG,n_convolution);
    phase = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));
    if all(phase == 0) || any(isnan(phase))
        error("")
    end

    % calculate real
    obsPAC = abs(mean(pwr.*exp(1i*phase)));

    % calculate nulls distribution
    for i=1:num_iter
        
        % Permutation method 1: select random time point
        random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
        timeshiftedpwr   = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
        permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
        
        % Permutation method 2: totally randomize power time series
        permutedPAC(2,i) = abs(mean(pwr(randperm(length(pwr))).*exp(1i*phase)));
        
        % Permutation method 3: FFT-based power time series randomization
        f = fft(pwr); % compute FFT
        A = abs(f);   % extract amplitudes
        zphs=cos(angle(f))+1i*sin(angle(f)); % extract phases
        powernew=real(ifft(A.*zphs(randperm(length(zphs))))); % recombine using randomized phases (note: use original phases to prove that this method reconstructs the original signal)
        powernew=powernew-min(powernew);
        
        permutedPAC(3,i) = abs(mean(powernew.*exp(1i*phase)));
    end

    % z-score real on nulls
    for i=1:3
        pacz(i,freq_idx) = (obsPAC-mean(permutedPAC(i,:)))/std(permutedPAC(i,:));
    end
end

figure;
for i=1:3
    switch i
        case 1
            % plot(EEG.times(1:end-1),timeshiftedpwr)
            title(['Null Method' num2str(i) ': Time-shifted'])
        case 2
            % plot(EEG.times(1:end-1),pwr(randperm(length(pwr))))
            title(['Null Method' num2str(i) ': randomized'])
        case 3
            plot(EEG.times(1:end-1),powernew)
            title(['Null Method' num2str(i) ': FFT-derived randomization'])
    end
    subplot(1,3,i)
    plot(frequencies,pacz(i,:))
    % title(["Null Method" num2str(i)])
    xlabel("Frequency (Hz)")
    ylabel("PACz")
end
sgtitle(sprintf("PACz by frequency (N_perm = %s) with windowing", num2str(num_iter)))
