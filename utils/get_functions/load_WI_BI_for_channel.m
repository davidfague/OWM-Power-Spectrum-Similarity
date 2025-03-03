function [EMS_WI, EMS_BI] = load_WI_BI_for_channel(params)
    [WI_channel_map] = load_BT_channel_map(params, 'WI');
    WI = WI_channel_map(params.chan_id);
    clear WI_channel_map
    [BI_channel_map] = load_BT_channel_map(params, 'BI');
    BI = BI_channel_map(params.chan_id);
    clear BI_channel_map
    [EMS_WI, EMS_BI] = pre_process_WI_BI(params, WI, BI);
end