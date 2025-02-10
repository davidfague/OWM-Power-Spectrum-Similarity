%% PSVs

% load
middle_fix = matfile('C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\OWM-Power-Spectrum-Similarity\processed_data\middle_fixation_baseline\201901\PSVs.mat');
end_fix = matfile('C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\OWM-Power-Spectrum-Similarity\processed_data\allpatients gammamod allregions allitem allenc baseline across trials\201901\PSVs.mat');

%%
row_to_compare = 10;

% % single trial
% PSV_whole_trials = {middle_fix.all_windowed_mean_PS_vectors(:,:,row_to_compare)', ...
%     end_fix.all_windowed_mean_PS_vectors(:,:,row_to_compare)'};

% % mean across trials
PSV_whole_trials = {mean(middle_fix.all_windowed_mean_PS_vectors(:,:,:),3)', ...
    mean(end_fix.all_windowed_mean_PS_vectors(:,:,:),3)'};

abs_max = max(max(PSV_whole_trials{1}(:), PSV_whole_trials{2}(:)));
abs_min = min(min(PSV_whole_trials{1}(:), PSV_whole_trials{2}(:)));

figure;
subplot(3,1,1)
imagesc(PSV_whole_trials{1})
title('middle fixation')
ylabel('frequency')
xlabel('time')
colorbar;
% clim(computeClims(PSV_whole_trials{2}))
clim([abs_min, abs_max])

subplot(3,1,2)
imagesc(PSV_whole_trials{2})
title("end fixation")
colorbar;
% clim(computeClims(PSV_whole_trials{2}))
clim([abs_min, abs_max])


subplot(3,1,3)
imagesc(PSV_whole_trials{1} - PSV_whole_trials{2})
title("middle - end")
colorbar;
% clim(computeClims(PSV_whole_trials{1} - PSV_whole_trials{2}))
% clim([abs_min, abs_max])


sgtitle("Baselining difference (PSVs) single trial")
%% EMS