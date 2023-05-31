function [dataI,bidx] = interpnan_NT(data)
    %linearly interpolate missing data points identified by NaNs
    % code given to me by Nathan Tardiff
    
    bidx = isnan(data); %bad (missing) data indicies
    gidx = ~bidx; %good data indicies

    %Interpolate
    dataI = data;
    %extrapolation necessary to handle possibility of NaNs at beginning/end
    %of data, extrapolating using median of data--this will work better if
    %preprocessing done on shorter segments (e.g., trials)
    %NO LONGER EXTRAP JUST CAUSES TROUBLE
    try
        %dataI(bidx) = interp1(find(gidx),dataI(gidx),find(bidx),'linear','extrap');
        %dataI(bidx) = interp1(find(gidx),dataI(gidx),find(bidx),'linear',nanmedian(dataI));
        dataI(bidx) = interp1(find(gidx),dataI(gidx),find(bidx));
    catch ME
        warning(ME.message);
    end
end