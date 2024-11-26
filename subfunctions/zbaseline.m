% % ADDED NOTE: Z-scoring is within channel - the mean, std is calculated with respect to individual channel across all trials.

function tfrz = zbaseline(data, baseline)

tfrz = data;
tfrz.powspctrm = zeros(size(tfrz.powspctrm));

ntimes = size(tfrz.powspctrm,4);

for e = 1:size(baseline.powspctrm,2)
    for f = 1:size(baseline.powspctrm,3)
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

end


% % Z-score TFR data on the baseline using statistical bootstrapping.
% %
% % E. L. Johnson, PhD
% % Copyright (c) 2017
% % UC Berkeley
% % eljohnson@berkeley.edu
% 
% function tfrz = zbaseline(data, baseline)
% 
% npermutes = 1000;
% 
% tfrz = data;
% tfrz.powspctrm = zeros(size(tfrz.powspctrm));
% 
% ntrials = size(tfrz.powspctrm,1);
% ntimes = size(tfrz.powspctrm,4);
% 
% for e = 1:size(baseline.powspctrm,2)
%     for f = 1:size(baseline.powspctrm,3)
%         base = squeeze(baseline.powspctrm(:,e,f,:));
%         base = base(:);
% 
%         basedist = zeros(1,npermutes);
%         for z = 1:npermutes
%             brand = randsample(base, ntrials);
%             basedist(z) = mean(brand); 
%         end
% 
%         bmean = repmat(mean(basedist), [ntrials 1 1 ntimes]);
%         bsd = repmat(std(basedist), [ntrials 1 1 ntimes]);
% 
%         tfrz.powspctrm(:,e,f,:) = (data.powspctrm(:,e,f,:)-bmean)./bsd;
% 
%         clear base basedist brand bmean bsd
%     end
% end
% 
% clearvars -except tfrz
% 
% end
