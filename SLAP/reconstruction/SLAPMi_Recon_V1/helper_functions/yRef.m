function yOut = yRef(yIn)
    %get an estimate of the F0 line intesities
    nCut = round(size(yIn,2)/4);
    [~, sortorder] = sort(nansum(yIn,1));
    yOut = nanmean(yIn(:, sortorder(10:10+nCut)), 2);
    yOut(sum(~isnan(yIn),2)<30) = nan;
    
    %replace nans with prediction from surrounding measurements
    yOut = max(0, inpaint_nans(yOut));
end
