function encoding_times = get_encoding_times_from_enc_ID(enc_ID, times, stim1_start, stim1_end, stim2_start, stim2_end, stim3_start, stim3_end)
% returns logical array
    if enc_ID == 1
        encoding_times = (times >= stim1_start & times < stim1_end);
    elseif enc_ID == 2
        encoding_times = (times >= stim2_start & times < stim2_end);
    elseif enc_ID == 3
        encoding_times = (times >= stim3_start & times < stim3_end);
    end
end