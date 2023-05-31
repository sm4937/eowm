function [betas] = run_pupil_GLM(alldata,sample_Hz)
%run_pupil_GLM Takes in all trials' timecourses of pupil size
%   alldata is ntrials x nsamples in size
%   sample_Hz is the eyelink sampling rate
%   betas is one per trial
%   constructs pupillary response function (PRF) by doing a trial-triggered
%   average, norming & smoothing it

PRF = nanmean(alldata,1);
unusable = isnan(PRF)|sum(isnan(alldata),1)>0;
PRF(unusable) = [];
PRF = PRF./max(abs([min(PRF) max(PRF)]));
% normalize like the gamma function would be (sum = 1, roughly)
% apply lowpass filter to PRF (smooth it)
PRF = lowpass(PRF,5,sample_Hz);
%specifies that (1) has been sampled at a rate of (3) hertz. (2) is the passband frequency of the filter in hertz.

% model = ones(size(alldata,1), length(PRF));
% % capture the number of correct trials and pupil samples in delay period
% %pres_model = zeros(size(pres_PRF));
% model = conv2(model,PRF,'valid');

modelInv = pinv(PRF);
alldata(:,unusable) = [];
betas = alldata * modelInv;


end

