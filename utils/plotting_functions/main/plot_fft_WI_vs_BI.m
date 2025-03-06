function fig = plot_fft_WI_vs_BI(WI,BI, in_dB, title_prefix, params)
    % WI and BI are of size (eTimes, mTimes, nTrialPairs)
    % Define colors (these are MATLAB's default "Lines" colors #1 and #2)
    WI_color = [0.0000, 0.4470, 0.7410];  % blue
    BI_color = [0.8500, 0.3250, 0.0980];  % reddish-orange

    srate = 1000/params.window_step_size; % convert from 1/ms to 1/s
    npoints = size(WI,2);
    f = (0:npoints-1)*(srate/npoints);       % Frequency vector
    f = f(1:npoints/2+1);           % Adjust frequency range

    % number of trial pairs
    num_WI = size(WI,3);
    num_BI = size(BI,3);

    % mean out encoding times, left with (1, maintenance times, trial pairs)
    WI = mean(WI,1);
    BI = mean(BI,1);

    % calculate num_frequencies
    num_frequencies = length(f);
    num_power_values = num_frequencies;

    WI_ffts = zeros(num_power_values, num_WI); %
    BI_ffts = zeros(num_power_values, num_BI);

    % compute WI_fftsand BI_ffts for each i in the 3rd dimensions of WI and
    % BI will have size (num_frequencies,
    % num_power_values=num_frequencies_values, and num_WI/num_BI)
    for i = 1:num_WI
        P1 = compute_fft(WI(:,:,i), srate, in_dB);
        WI_ffts(:,i) = P1;
    end
    for i = 1:num_BI
        P1 = compute_fft(BI(:,:,i), srate, in_dB);
        BI_ffts(:,i) = P1;
    end


    %% subplot(311) plot the mean with shaded SEM for WI_ffts and BI_ffts
    mean_WI = squeeze(mean(WI_ffts, 2));
    mean_BI = squeeze(mean(BI_ffts, 2));
    WI_err = squeeze(std(WI_ffts, 0, 2) / sqrt(size(WI_ffts, 2)));     % Standard deviation (0 flag for normalization by N-1) %  / sqrt(n) to convert to standard error of the mean
    BI_err = squeeze(std(BI_ffts, 0, 2) / sqrt(size(BI_ffts, 2)));

    % Calculate the upper and lower bounds of the error region
    mean_WI_upper = mean_WI + WI_err;
    mean_WI_lower = mean_WI - WI_err;
    mean_BI_upper = mean_BI + BI_err;
    mean_BI_lower = mean_BI - BI_err;

    % plot
    fig = figure;
    subplot(1,1,1); title('mean'); hold on; % means
    legend("AutoUpdate","on")
    fill([f, fliplr(f)], [mean_WI_upper', fliplr(mean_WI_lower')], ...
        WI_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, ... % Make the fill semi-transparent
        'DisplayName', 'WI SEM');
    plot(f, mean_WI, 'DisplayName', 'Within-Images')
    % --- Between-Images shaded area & line ---
    fill([f, fliplr(f)], ...
         [mean_BI_upper', fliplr(mean_BI_lower')], ...
         BI_color, ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none', ...
         'DisplayName', 'BI SEM');
    plot(f, mean_BI, 'DisplayName', 'Between-Images')
    % % subplot(312) plot each of the WI_ffts in the third dimension
    % subplot(3,1,2); % individual trials
    % plot(f, WI_ffts)
    % title('WI trials')
    % 
    % % subplot(313) plot each of the BI_ffts in the third dimension
    % subplot(3,1,3); % individual trials
    % plot(f, BI_ffts)
    % title('BI trials')

    title_str = sprintf("%s p%s chan%s image%s enc%s\n %s", ...
            params.type, num2str(params.patient_id), ... 
            num2str(params.chan_id), num2str(params.image_id), ...
            num2str(params.enc_id), params.anat);
    sgtitle(title_str);

    xlabel('Frequency (Hz)')
    if in_dB
        ylabel('Magnitude (dB)')
    else
        ylabel('Magnitude')
    end
    title(sprintf('%s Amplitude Spectrum', title_prefix))

end


function P1 = compute_fft(signal_to_fft, srate, in_dB)
    fft_signal = fft(signal_to_fft);
    npoints = length(signal_to_fft);
    % f = (0:npoints-1)*(srate/npoints);       % Frequency vector
    P2 = abs(fft_signal/npoints);            % Normalize FFT output
    P1 = P2(1:npoints/2+1);         % Take only the positive half (since FFT is symmetric)
    % f = f(1:npoints/2+1);           % Adjust frequency range
    if in_dB
        P1 = 20*log10(abs(P1));
    end
end