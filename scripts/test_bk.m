%% EMS trace analysis with xcorr


load("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc correct\201907\encoding_similarity.mat")
load("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc correct\201907\PSVs.mat", 'label_table')
% load("D:\Power Spectrum Similarity\Z_Raw Data Storage\D_OWM_t_bipolar_201907.mat")
load("D:\Power Spectrum Similarity\Z_Raw Data Storage\OWM_trialinfo_201907.mat", 'C')
rows_without_nan_by_enc = ~squeeze(any(any(isnan(all3_ES_matrix), 2), 3));
clearvars -except all3_ES_matrix label_table rows_without_nan_by_enc C
%% select channel, correctness, plot ems avg across trials for each encoding

chan=66;
% targID=3;
correctness=1;
Mtime=[250:600]; % delay time
for targID=1:9
    figure;
    sgtitle(C{1,targID});
    all3_data = [];
    trials =[];
    for i=1:3
        subplot(1,3,i)
        encodingperiod=i;

        clear targtable data
        targtable=label_table(label_table.channel_ID==chan & label_table.encID_to_imageID(:,encodingperiod)==targID & label_table.encoding_correctness(:,encodingperiod)==correctness ,:);
        % targtable=strcmp(targtable.anatomical_label, "Hippocampus"
        data=all3_ES_matrix(label_table.channel_ID==chan & label_table.encID_to_imageID(:,encodingperiod)==targID & label_table.encoding_correctness(:,encodingperiod)==correctness,:,:,:);

        imagesc(squeeze(mean(data(:,:,Mtime,i),1, 'omitnan')), [0 .5])

        EMS_to_plot = squeeze(data(:,:,Mtime,i));
        all3_data = cat(1, all3_data, EMS_to_plot);

    end
    % subplot(1,4,4)
    % imagesc(squeeze(mean(data(:,:,Mtime,:),1, 'omitnan'), 4), [0 .5])
    % subplot(1,4,4);
    % imagesc(flipud(squeeze(mean(all3_data(:,:,:),1, 'omitnan'))), [0 .5])

end
% subplot(1,4,4);
% imagesc(flipud(squeeze(mean(all3_data(:,:,:),1, 'omitnan'))), [0 .5])
        %% checking recomputing with -0.5 to 0 baseline
        % targtable = label_table(label_table.channel_ID==chan & label_table.encID_to_imageID(:,encodingperiod)==targID & label_table.encoding_correctness(:,encodingperiod)==correctness,:);
        % for row_idx = 1:size(targtable,1)
        %     channel_id = targtable(row_idx,2);
        %     trial_id = targtable(row_idx,end);
        % 
        %     compute_Zpower(D_OWM_t(channel_id,:,trial_id))
        
        %%
% figure;
% imagesc(squeeze(mean(mean(data,1),4) ),[0 .5])

%% trial by trial

chan=1;
targID=1;
correctness=1;
encodingperiod=3;


clear data targtabletr
targtabletr=label_table(label_table.channel_ID==chan & label_table.encID_to_imageID(:,encodingperiod)==targID & label_table.encoding_correctness(:,encodingperiod)==correctness ,:);
data=all3_ES_matrix(label_table.channel_ID==chan & label_table.encID_to_imageID(:,encodingperiod)==targID & label_table.encoding_correctness(:,encodingperiod)==correctness,:,:,:);

figure;
sgtitle(C{1,targID});
for i=11:20 %size(targtabletr,1)

    subplot(10,1,i-10)



    imagesc(squeeze(mean(data(i,:,:,encodingperiod),1, 'omitnan')), [0 .5])


end

%% look at ES avg across enc over time
% clearvars -except all3_ES_matrix label_table
chan=50;

correctness=1;
figure;
% sgtitle(C{1,targID});

% rowID=find(label_table.channel_ID==chan  & sum(label_table.encoding_correctness,2)==3);
rowID=find(label_table.channel_ID==chan  & sum(label_table.encoding_correctness,2)==3 & sum(rows_without_nan_by_enc, 2)==3);

targtable = label_table(rowID,:);

for i=1:5 % trial
subplot(5,1,i)

    for encodingperiod=1:3
        plot(squeeze(mean(all3_ES_matrix(rowID(i),:,:,encodingperiod),2, 'omitnan'))), hold on
        legend()
    end
end
%% check same ES different trials
label_table(rowID(3:5),:) % something definitely wrong bc same means

mf = matfile("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc\201907.mat");
mf.all_windowed_mean_PS_vectors(1,1, rowID(3:5)) % they are the same, unfortunately

% check channel 10, trial 19:21 201907
D_OWM_t = matfile("D:\Power Spectrum Similarity\Z_Raw Data Storage\D_OWM_t_bipolar_201907.mat");
D_OWM_t.D_OWM_t(10,8900:8905,19:21)

sum(isnan(D_OWM_t.D_OWM_t(:,:,:)), [1 2 3])

% check onedrive verson
D_OWM_t = matfile("C:\Users\drfrbc\Downloads\D_OWM_t_bipolar.mat");
D_OWM_t.D_OWM_t(10,8900:8905,19:21)

% checking 201908
D_OWM_t = matfile("C:\Users\drfrbc\Downloads\D_OWM_t_bipolar (1).mat");
D_OWM_t.D_OWM_t(10,8900:8905,19:21)

load("C:\Users\drfrbc\Downloads\D_OWM_t_bipolar (1).mat", 'D_OWM_t');
% for channelID = 1:90
%     for trial1 = 1:180
%         for trial2 = 1:180
%             if D_OWM_t(channelID,2,trial1)==D_OWM_t(channelID,2, trial2)
%                 error()
%             end
%         end
%     end
% end

% not in 201908, 201906, 201905, 201903, 201902, 201901, 201910, 201915

load("C:\Users\drfrbc\Downloads\D_OWM_t_bipolar (5).mat", 'D_OWM_t');

% Extract the data for all trials of interest
data = squeeze(D_OWM_t(:, 5, :)); % Size will be (channels x trials)

% Check for duplicate values in each channel across trials
for channelID = 1:size(data, 1)
    if numel(unique(data(channelID, :))) < size(data, 2)
        fprintf('Duplicate values found in channel %d \n', channelID);
    end
end




%% single value for each trial
peaks =[];
delays =[];
corr12=[];
lag12=[];
Mtime=[250:600]; % delay time
% Mtime=[1:75]; fixation time

for i =  1:length(rowID)
    trace1 = squeeze(mean(all3_ES_matrix(rowID(i),:,Mtime,1),2));
    trace2 = squeeze(mean(all3_ES_matrix(rowID(i),:,Mtime,3),2));

    trace3 = squeeze(mean(all3_ES_matrix(rowID(i),:,:,3),2));


    [corr12(i,:), lag12(i,:)] = xcorr(trace1, trace2, 'coeff'); % Cross-correlation between trace1 and trace2
    % [peak, idx] = max(corr12); % Use max to get the highest peak
    % peaks(i) = peak;
    % delays(i) = lag12(idx);

end
% figure;
% plot(peaks,delays)

figure; plot(mean(lag12)',mean(corr12)')

figure; plot((lag12)',(corr12)')

%% look at cross correlation

% trace1 = squeeze(mean(all3_ES_matrix(rowID(i),:,:,1),2));
% trace2 = squeeze(mean(all3_ES_matrix(rowID(i),:,:,2),2));
trace3 = squeeze(mean(all3_ES_matrix(rowID(i),:,:,3),2));

trace1 = squeeze(mean(all3_ES_matrix(rowID(2),:,:,1),2));
trace2 = squeeze(mean(all3_ES_matrix(rowID(2),:,:,3),2));

% %replace inf with max
% trace1(isinf(trace1)) = max(trace1(~isinf(trace1)));
% trace2(isinf(trace2)) = max(trace2(~isinf(trace2)));

% get rid of everything before maintenance
lastinf = find(isinf(trace3), 1, 'last');
trace1 = trace1(lastinf:end);
trace2 = trace2(lastinf:end);



[corr12, lag12] = xcorr(trace1, trace2, 'coeff'); % Cross-correlation between trace1 and trace2

% Plot the cross-correlations
figure;

% Plot trace1 vs trace2 cross-correlation
subplot(3,1,1);
plot(lag12, corr12);
title('Cross-Correlation: Trace 1 vs Trace 2');
xlabel('Lag');
ylabel('Correlation Coefficient');


%% look at cross correlation over time

% Parameters for binning
bin_size = 10; % Bin size (10 points) (each point is a window ID)
n_points = length(trace1);
n_bins = floor(n_points / bin_size); % Total number of bins

% Initialize storage for peaks and lags
peak_values = zeros(n_bins, 1);
peak_lags = zeros(n_bins, 1);

lags_bins = (-bin_size + 1):(bin_size - 1); % Lag values are constant for all bins

% Loop over bins
for b = 1:n_bins
    % Define bin range
    bin_start = (b - 1) * bin_size + 1;
    bin_end = bin_start + bin_size - 1;

    % Extract data for the current bin
    bin_trace1 = trace1(bin_start:bin_end);
    bin_trace2 = trace2(bin_start:bin_end);
    % bin_trace3 = trace3(bin_start:bin_end);

    % Compute cross-correlation for the current bin
    [corr_bin, lag_bin] = xcorr(bin_trace1, bin_trace2, 'coeff');

    % Find the peak value and its lag
    [peak, idx] = max(corr_bin); % Use max to get the highest peak
    peak_values(b) = peak;
    peak_lags(b) = lag_bin(idx);

end

% Visualize the peak values and delays over time
time_bins = (1:n_bins) * bin_size; % Time for each bin
figure;

% Plot peak values
subplot(2, 1, 1);
plot(time_bins, peak_values, 'o-', 'LineWidth', 1.5);
title('Peak Cross-Correlation Over Time');
xlabel('Time (points)');
ylabel('Peak Value');
grid on;

% Plot peak lags
subplot(2, 1, 2);
plot(time_bins, peak_lags, 'o-', 'LineWidth', 1.5);
title('Lag of Peak Cross-Correlation Over Time');
xlabel('Time (points)');
ylabel('Lag (points)');
grid on;

%% check nans
% Example matrix size: [9460, 41, 640, 3]
% Assume all3_ES_matrix is loaded or already in the workspace

% Logical array indicating where NaNs are present
nan_locations = isnan(all3_ES_matrix);

% Collapse the second and third dimensions by checking for any NaNs
% across these dimensions
nan_summary = squeeze(any(any(nan_locations, 2), 3)); % Size: [9460, 3]

% Visualize the NaN summary
figure;
imagesc(nan_summary'); % Transpose for better visualization
colormap('hot'); % Use a heatmap color scheme
colorbar;
xlabel('First Dimension (Index)');
ylabel('Fourth Dimension (Index)');
title('NaN Locations Across First and Fourth Dimensions');