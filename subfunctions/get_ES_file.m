function ES_file = get_ES_file(PS_file)
    ES_file = matfile(fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('all3_ES.mat'))));
    if ~exist(ES_file.Properties.Source, 'file')
        ES_file = matfile(fullfile(strrep(PS_file.Properties.Source, '.mat', 'all3Enc_ES.mat')));
    end
    
    if ~exist(ES_file.Properties.Source, 'file')
        error("ES file does not exist: %s", ES_file.Properties.Source)
    end
end