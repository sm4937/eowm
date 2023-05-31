function [betas,residuals] = deconvolve_pupil(Y,model)
% deconvolve_pupil takes in pupil size data, eye tracking data, and
% frees the pupil timecourse from confounds
%   data is already cleaned, band-pass filtered pupil size
%   residuals are the output, truly squeaky clean pupil size
%   timecourses

    % betas = pinv(model)*Y;
    [betas,BINV,~,~,stats] = regress(Y,model);
    % get effect of intercept, difficulty condition, X, Y eye position
    % & baseline

    % predict data with betas (but not difficulty betas because then we
    % can't see condition differences)
    model(:,2) = []; betas_sim = betas; betas_sim(2) = [];
    Y_hat = model*betas_sim;
    % reconstruct pupil time course w/out X/Y or baseline information

    residuals = Y - Y_hat;
    % this is your new delay timecourse
    residuals = zscore(residuals);

end

