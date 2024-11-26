# OWM-Power-Spectrum-Similarity
Analyzes iEEG human brain signals recorded during an Object Working Memory task by calculating similarities between time-binned frequency-power vectors.

Abreviated Terminology:
Power Spectrum Similarity (PSS)
Power Spectrum Vectors (PSVs)

Data:
  (1) Raw data is intracranial electroencephalgram signals recorded from ~10 patients, ~70-90 channels each, ~100-300 trials each. 0-1000ms fixation, 1000-1500ms encoding1, 1500-2000ms encoding2, 2000-2500ms encoding3, 2500-6500ms maintenance, 6500-9000ms recall.
  (2) PSVs include frequencies 1:40 hz. PSV time windows are 100 ms wide, with 10 ms step. (methods)

Objectives:
  (1) Compute PSVs for each patient, channel, trial; baseline and z-score to fixation.
  (2) Compute the within-trial similarity for every patient, channel, trial, and all 3 objects encoded in the trial.
  (3) Compute the between-trial similarity between items.
  (4) Compute the between-trial similarity within items.

Workflow:
  (1)
