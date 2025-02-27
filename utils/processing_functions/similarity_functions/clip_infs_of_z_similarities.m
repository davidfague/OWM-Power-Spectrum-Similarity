function [WI, BI] = clip_infs_of_z_similarities(WI, BI)
    if isempty(WI(isinf(WI))) && isempty(BI(isinf(BI)))
        return % no infs to clip
    else
        best_z = find_best_z_clipping(WI, BI);
        bounds = [-best_z, best_z];
        
        WI(WI>bounds(2)) = bounds(2);
        WI(WI<bounds(1)) = bounds(1);
    
        BI(BI>bounds(2)) = bounds(2);
        BI(BI<bounds(1)) = bounds(1);
    end
end