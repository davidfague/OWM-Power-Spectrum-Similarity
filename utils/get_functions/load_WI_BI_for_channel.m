function [EMS_WI, EMS_BI] = load_WI_BI_for_channel(params)
    [WI_channel_map] = load_BT_channel_map(params, 'WI');
    try
        WI = WI_channel_map(params.chan_id);
    catch ME
        warning('chan %s is not in WI_channel_map', params.chan_id)
        EMS_WI = [];
        EMS_BI = [];
        return
    end
    clear WI_channel_map
    [BI_channel_map] = load_BT_channel_map(params, 'BI');
    try
        BI = BI_channel_map(params.chan_id);
    catch ME
        warning('chan %s is not in BI_channel_map', params.chan_id)
        EMS_WI = [];
        EMS_BI = [];
        return
    end
    clear BI_channel_map
    [EMS_WI, EMS_BI] = pre_process_WI_BI(params, WI, BI);
end