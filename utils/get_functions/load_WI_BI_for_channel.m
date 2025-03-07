function [EMS_WI, EMS_BI] = load_WI_BI_for_channel(params)
    [WI_channel_map] = load_BT_channel_map(params, 'WI');
    if isempty(WI_channel_map)
        EMS_WI = [];
        EMS_BI = [];
        return
    else
        WI = WI_channel_map(params.chan_id);
    end
    clear WI_channel_map
    [BI_channel_map] = load_BT_channel_map(params, 'BI');
    if isempty(BI_channel_map)
        EMS_WI = [];
        EMS_BI = [];
        return
    else
        BI = BI_channel_map(params.chan_id);
    end
    clear BI_channel_map
    [EMS_WI, EMS_BI] = pre_process_WI_BI(params, WI, BI);
end