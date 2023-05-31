function [filtered_and_cleaned_data] = clean_pupil_data(data,sample_Hz)
%clean_pupil_data Takes in pupil data, of any length, and returns bandpass
%filtered & cleaned (interpolated) version of it
% Assumes one trial (one row/column) of data at a time.
% data is the input data
% sample_Hz is the sampling rate of the original data
% filtered_and_cleaned_data is the cleaned version

    cleaned_data = interpnan_NT(data); 
    % linearly interpolate data

    % filter
    % determine if any NaNs left over at edges after interpolation (no longer
    % extrapolating due to data distortion)
    edgeNaNs = isnan(cleaned_data);

    % apply lowpass butterworth filter
    options.filter_cutoff = 8;
    options.filter_order = 3;

    norm_cutoff = options.filter_cutoff./(sample_Hz./2); %normalized cutoff freq
    [b,a] = butter(options.filter_order,norm_cutoff);
    try
        filtered_and_cleaned_data = nan(size(cleaned_data));
        filtered_and_cleaned_data(~edgeNaNs) = filtfilt(b,a,cleaned_data(~edgeNaNs));
    catch ME
        %handles cases where there are too few valid samples to actually
        %filter, which will presumably be thrown away at QA anyway
        warning(ME.message);
        filtered_and_cleaned_data = cleaned_data;
    end
    %another sanity check
    assert(isequal(edgeNaNs,isnan(filtered_and_cleaned_data)),'Filtering induced NaNs!')
    
end

