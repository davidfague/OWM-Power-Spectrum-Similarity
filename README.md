# OWM-Power-Spectrum-Similarity
Analyzes human iEEG brain signals recorded during an Object Working Memory task by calculating similarities between time-binned frequency-power vectors.

Abreviated Terminology:
Power Spectrum Similarity (PSS)
Power Spectrum Vectors (PSVs)
Within-Image (WI)
Between-Image (BI)

Data:
  (1) Raw data is pre-processed intracranial electroencephalgram (iEEG) signals recorded from ~20 patients, ~70-150 channels each, ~100-300 trials each. 0-1000ms fixation, 1000-1500ms encoding1, 1500-2000ms encoding2, 2000-2500ms encoding3, 2500-6500ms maintenance, 6500-9000ms recall.
  (2) Data will not be made available until after publication
  (3) PSVs include frequencies 1:40 hz (adjustable). PSV time windows are 100 ms wide, with 10 ms step (adjustable). (see methods)
  (4) similarities are calculated using Spearman-rank correlation on PSVs (see methods)

Objectives:
  (1) Compute PSVs for each patient, channel, trial; baseline and z-score to fixation.
  (2) Compute the within-trial similarity for every patient, channel, trial, and all 3 objects encoded in the trial.
  (3) Compute the between-trial similarity between items.
  (4) Compute the between-trial similarity within items.

Workflow:
  Scripts:
    (1) get_parameters.m:
              Overview: function containing default parameters for quick initialization/sharing/adjustment across scripts. 
                Contains: pathing and processing hyperparameters
    (2) zpower_and_psvs.m:
              Overview: converts pre-processed data into PSVs of baseline-corrected trial-induced, z-scored power across specified frequencies.
                Slices and manages pre-processed data, calls compute_Zpower on the data slice, and saves a large matrix of PSVs (all_windowed_mean_PS_vectors) and info table (label_table) of lengths = (channels*trials)
              Calls helper functions: 
                compute_Zpower which calls:
                    BOSC_tf_power() for computing power using Morelet Wavelet approximation
                    zbaseline() for z-scoring power to 'fixation period'
    (3) compute_es_between_trials.m
              Overview: converts PSVs into similarities between trial pairs. 
                Notes: computes similarities between 'encoding period' and either 'maintenance period' or 'encoding period' between trials.
                       'between trials' has trial pairing types: 'within image' (WI) or 'between image' (BI)
                       same-trial 'within image' pairs are left in until post-processing visualization
              Calls processing function:
                  compute_all_between_trial_similarities: 
                      finds the data slices
                      for each channel
                        computes a large matrix ("BT_ES") of size = (time_windows1, time_windows2, ctrl_trials, test_trials) and values: rho similarity
                        computes info tables ("chan_test_table", "chan_ctrl_table") containing the info for the 'ctrl_trials, test_trials' slices
              helper functions:
                  compute_BT_similarity_matrix:
                      computes the similarity matrix (time1, time2) for data slice (trial1, trial2)

    (4) compute_es_within_trials.m
          Overview: same as "compute_es_between_trials.m" except processes pairs each row of the all_windowed_mean_PS_vectors with itself.
          computes 'encoding vs whole-trial' instead of 'encoding vs maintenance or encoding'

    (5) main_plotting.m
          Overview: visualizes the results
          (a) plot_PSVs
