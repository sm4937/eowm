%% A quick script to run power analyses on TAFKAP outputs
%Sarah Master, EOWM 2021
clear all

subjs = [5 6];%[4 5 6];
ROI = 'IPS0';
estimates = []; uncertainty = []; stimvalues = []; conds = [];
longest_length = 240;

%Look into pilot subject (number 4)
for sii = 1:length(subjs)
    
    subjnum = subjs(sii);
    load(['TAFKAP_subj' num2str(subjnum) '_late_delay.mat'])
    load(['task_subj' num2str(subjnum) '.mat'])
    
    temp = NaN(longest_length,1);
    eval(['temp_estimates = TAFKAP_output.' ROI '.est;'])
    temp(1:length(temp_estimates),:) = temp_estimates;
    estimates = [estimates temp];
    
    temp = NaN(longest_length,1);
    eval(['temp_unc = TAFKAP_output.' ROI '.unc;'])
    temp(1:length(temp_unc),:) = temp_unc;
    uncertainty = [uncertainty temp];
    
    temp = NaN(longest_length,1);
    eval(['temp_stims = task_subj' num2str(subjnum) '.stimval(~isnan(task_subj' num2str(subjnum) '.stimval))'])
    temp(1:length(temp_stims),:) = temp_stims;
    stimvalues = [stimvalues temp];
    
    temp = NaN(longest_length,1);
    eval(['temp_conds = task_subj' num2str(subjnum) '.cond(~isnan(task_subj' num2str(subjnum) '.stimval))'])
    temp(1:length(temp_conds),:) = temp_conds;
    conds = [conds temp];
end

error = get_angular_distance(estimates.*2, stimvalues);

for sii = 1:length(subjs)
    
    uncmeans(sii,:) = [nanmean(uncertainty(conds(:,sii)==1,sii)) nanmean(uncertainty(conds(:,sii)==2,sii))];
    uncstds(sii,:) = [nanstd(uncertainty(conds(:,sii)==1,sii)) nanstd(uncertainty(conds(:,sii)==2,sii))];
    
    errormeans(sii,:) = [nanmean(error(conds(:,sii)==1,sii)) nanmean(error(conds(:,sii)==2,sii))];
    errorstds(sii,:) = [nanstd(error(conds(:,sii)==1,sii)) nanstd(error(conds(:,sii)==2,sii))];
    
end


figure
subplot(1,2,1)
errorbar(nanmean(errormeans),nanmean(errorstds)./sqrt(length(subjs)),'ok','LineWidth',2)
title('Decoding error across conditions')
ylabel('Error')

subplot(1,2,2)
errorbar(nanmean(uncmeans),nanmean(uncstds)./sqrt(length(subjs)),'ok','LineWidth',2)
title('Decoding uncertainty across conditions')
ylabel('Uncertainty')
fig = gcf; fig.Color = 'w';

% I think this needs to be done w more than 1 subject
% So take the difference between condition means for each subject
% Take the std of subjects' difference scores
% That is the std of the null & the alternative hypothesis distributions
% The null mean is 0, no difference in conditions
% The alternative hypothesis mean is the real sample mean difference
% So mean(errordiffs) with each subject providing one score

N = sampsizepwr('t',[0 std(errordiffs)], mean(errordiffs), 0.95)
% enter test type (1-sample t-test (paired)), then mean & std of null
% hypothesis (no difference), and mean of alt hypothesis (yes difference),
% then desired power

% Another way of doing it, which makes less sense I think?
N = sampsizepwr('t',[uncmeans(1) uncstds(1)], uncmeans(2),0.99)

