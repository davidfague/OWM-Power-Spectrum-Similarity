% ADDED NOTE: Z-scoring is within channel - the mean, std is calculated with respect to individual channel across all trials.
% powspctrm is (channels, trials, frequencies, times)
function tfrz = zbaseline(data, baseline, baseline_across_trials)

if ~baseline_across_trials

    tfrz = data;
    tfrz.powspctrm = zeros(size(tfrz.powspctrm));
    
    ntimes = size(tfrz.powspctrm,4);
    
    for e = 1:size(baseline.powspctrm,2) % trials
        for f = 1:size(baseline.powspctrm,3) % frequencies
            base = squeeze(baseline.powspctrm(:,e,f,:));
            base = base(:);
    
            % Directly calculate the mean and std from the baseline for ntrials = 1
            bmean = repmat(mean(base), [1 1 1 ntimes]);
            bsd = repmat(std(base), [1 1 1 ntimes]);
    
            % Perform z-scoring
            tfrz.powspctrm(:,e,f,:) = (data.powspctrm(:,e,f,:) - bmean) ./ bsd;
    
            clear base bmean bsd
        end
    end
    
    clearvars -except tfrz

% end
else
% original, across trials method:
% Z-score TFR data on the baseline using statistical bootstrapping.
%
% E. L. Johnson, PhD
% Copyright (c) 2017
% UC Berkeley
% eljohnson@berkeley.edu

% function tfrz = zbaseline(data, baseline)

npermutes = 1000;

tfrz = data;
tfrz.powspctrm = zeros(size(tfrz.powspctrm));

ntrials = size(tfrz.powspctrm,1);
ntimes = size(tfrz.powspctrm,4);

for e = 1:size(baseline.powspctrm,2) % channels
    for f = 1:size(baseline.powspctrm,3) %frequencies
        base = squeeze(baseline.powspctrm(:,e,f,:)); % the baseline's time-power data from all trials for this (channel,frequency)
        base = base(:);

        basedist = zeros(1,npermutes);
        for z = 1:npermutes
            brand = randsample(base, ntrials); % pull ntrials random values from all 'base' values
            basedist(z) = mean(brand); % mean
        end

        bmean = repmat(mean(basedist), [ntrials 1 1 ntimes]); % compute mean of basedist
        bsd = repmat(std(basedist), [ntrials 1 1 ntimes]); % compute std of basedist

        tfrz.powspctrm(:,e,f,:) = (data.powspctrm(:,e,f,:)-bmean)./bsd; % for all values subtract mean and divide by standard deviation

        clear base basedist brand bmean bsd
    end
end

clearvars -except tfrz

end
