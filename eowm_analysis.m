%% Data analysis script for EOWM
% Sarah Master
% 01/2021
clear all
% first-pass at analyzing effort of WM task data
files = dir('data/'); files = char(files.name);
data_folders = files(contains(string(files),'subj'),:);

addpath(genpath('../superbar_all_files'))

subjs = [4:14 16]; 

data = [];
for s = 1:length(subjs)
    subj = subjs(s);
    data = [data; makeEOWMdatatable(subj)];
end
n = length(data.subj); condlabels = {'Easy','Hard'}; 
condcolors = [190 0 110; 0 110 190]./255; condcolors = flip(condcolors,1);
% Hard is pink, easy is blue, keep it that way

disp(['Mean delta for easy trials: ' num2str(mean(data.deltas(:,1)))])
disp(['Mean delta for hard trials: ' num2str(mean(data.deltas(:,2)))])

disp(['Mean RT for easy trials: ' num2str(mean(data.cond_rt(:,1)))])
disp(['Mean RT for hard trials: ' num2str(mean(data.cond_rt(:,2)))])

disp(['Between ' num2str(100*min(data.pct_excluded)) '% and ' num2str(100*max(data.pct_excluded)) '% trials excluded per participant.'])
disp(['Mean(std): ' num2str(mean(data.pct_excluded)) '(' num2str(std(data.pct_excluded)) ')']);

pupil.subj = subjs';
%pupil.delay_betas = data.delay_betas;
%pupil.stim_betas = data.stim_betas;
pupil.trial_indices = data.pupil_trial_indices;
save('data/pupil.mat','pupil')

for subj = 1:length(subjs)
    
    delay_pupils = data.cleaned_pupilsize_all_delays{subj};
    hard_timecourses_delay(subj,:) = nanmean(delay_pupils(delay_pupils(:,1)==2,2:end));
    easy_timecourses_delay(subj,:) = nanmean(delay_pupils(delay_pupils(:,1)==1,2:end));
    
    stim_pupils = data.cleaned_pupilsize_all_stims{subj};
    hard_timecourses_stim(subj,:) = nanmean(stim_pupils(stim_pupils(:,1)==2,2:end));
    easy_timecourses_stim(subj,:) = nanmean(stim_pupils(stim_pupils(:,1)==1,2:end));
    
    wholetrial_pupils = data.cleaned_pupilsize_by_trial{subj};
    hard_timecourses_trial(subj,:) = nanmean(wholetrial_pupils(wholetrial_pupils(:,1)==2,2:end));
    easy_timecourses_trial(subj,:) = nanmean(wholetrial_pupils(wholetrial_pupils(:,1)==1,2:end));
    
end

data.mean_delay_pupil_timecourses_hard = hard_timecourses_delay;
data.mean_delay_pupil_timecourses_easy = easy_timecourses_delay;
data.mean_stim_pupil_timecourses_hard = hard_timecourses_stim;
data.mean_stim_pupil_timecourses_easy = easy_timecourses_stim;


%% Make plots appropriate for VSS Poster 2022

% Condition accuracy, withOUT t-tests
% Not averaged across subjects, but showing all subjects' behavior

figure
subplot(3,1,1)
plot(1:2,data.cond_accuracy,'--','LineWidth',0.5,'Color','k')
hold on
errorbar(1, mean(data.cond_accuracy(:,1)),nanstd(data.cond_accuracy(:,1))./sqrt(n),'o','LineWidth',3,'Color',condcolors(1,:),'DisplayName','Easy trials')
errorbar(2, mean(data.cond_accuracy(:,2)),nanstd(data.cond_accuracy(:,2))./sqrt(n),'o','LineWidth',3,'Color',condcolors(2,:),'DisplayName','Hard trials')
[h,pval] = ttest(data.cond_accuracy(:,1),data.cond_accuracy(:,2));
if pval < 0.05
    scatter(1.5,0.95,'k*','LineWidth',1.5)
end
xticklabels(condlabels)
xticks([1:2]); xlim([0.75 2.25]); 
plot(xlim,[0.7 0.7],'-.','LineWidth',1,'Color',condcolors(2,:),'DisplayName','Target accuracy for hard trials')
plot(xlim,[0.9 0.9],'-.','LineWidth',1,'Color',condcolors(1,:),'DisplayName','Target accuracy for easy trials')
[fig,ax] = clean_fig();
ylabel('% correct');

subplot(3,1,2)
plot(1:2,data.cond_rt,'--','LineWidth',0.5,'color','k')
hold on
errorbar(1, mean(data.cond_rt(:,1)),nanstd(data.cond_rt(:,1))./sqrt(n),'o','LineWidth',3,'Color',condcolors(1,:),'DisplayName','Easy trials')
errorbar(2, mean(data.cond_rt(:,2)),nanstd(data.cond_rt(:,2))./sqrt(n),'o','LineWidth',3,'Color',condcolors(2,:),'DisplayName','Hard trials')
[h,pval] = ttest(data.cond_rt(:,1),data.cond_rt(:,2));
if pval < 0.05
    scatter(1.5,0.59,'k*','LineWidth',1.5)
end
xticklabels(condlabels)
xticks([1:2]); xlim([0.75 2.25]); 
[fig,ax] = clean_fig();
ylabel('Mean RT (seconds)')

subplot(3,1,3)
fulltc = [nanmean(easy_timecourses_trial); nanmean(hard_timecourses_trial)];
full_SEM = [nanstd(easy_timecourses_trial,1);nanstd(hard_timecourses_trial,1)]./sqrt(n);
xs = [1:length(easy_timecourses_trial)].*(2/1000);

toclean = sum(isnan(fulltc),1)>0 & sum(isnan(full_SEM),1)>0;
fulltc(:,toclean) = []; full_SEM(:,toclean) = [];
xs(:,toclean) = [];

btwn_fill = [fulltc(2,:) + full_SEM(2,:), fliplr(fulltc(2,:))-fliplr(full_SEM(2,:))];
fill_xs = [xs, fliplr(xs)];
fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials');
hold on
btwn_fill = [fulltc(1,:) + full_SEM(1,:), fliplr(fulltc(1,:))-fliplr(full_SEM(1,:))];
fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials');
plot(xs,fulltc(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
%xlim([0 12])
legend({'Hard','Easy'})
xlabel('Time (seconds)')
ylabel('Z-scored pupil size')
[fig,ax] = clean_fig();

[h,p] = ttest2(nanmean(easy_timecourses_delay,2),nanmean(hard_timecourses_delay,2));
% is there a significant difference over whole delay period?
% no. p = 0.6
% nothing in early or late delay period, either
[h,p] = ttest2(nanmean(easy_timecourses_delay(:,1:3000),2),nanmean(hard_timecourses_delay(:,1:3000),2));
[h,p] = ttest2(nanmean(easy_timecourses_delay(:,3001:6000),2),nanmean(hard_timecourses_delay(:,3001:6000),2));

[h,p] = ttest2(nanmean(easy_timecourses_stim,2),nanmean(hard_timecourses_stim,2));
% is there a significant difference over stimulus presentation period?
% no. p = 0.4


%% Let's do a new pupil size analysis - cleaner, more principled

figure
title('Effect of hard trial on pupil size')
errorbar([nanmean(data.delay_betas(:,2)) nanmean(data.stim_betas(:,2))], ...
    [nanstd(data.delay_betas(:,2)) nanstd(data.stim_betas(:,2))]./sqrt(n), ...
    '*k','LineWidth',2.5)
xticklabels({'Delay period','Stimulus presentation period'})
xticks([1 2])
xlim([0.75 2.25])
ylabel('Regression \beta')
clean_fig();

[h,p] = ttest(data.delay_betas(:,2));

%what's up with the 0 beta weights for x-y position? 
% NaNs?


% In the style of Wiehler et al., 2022
% % regression-based, cleaner, within-subject-type analysis
% figure
% subplot(2,1,1)
% bar([nanmean(data.easy_delay_beta), nanmean(data.hard_delay_beta)],'w')
% hold on
% scatter(ones(n,1),data.easy_delay_beta,[],condcolors(1,:),'Filled')
% scatter(2.*ones(n,1),data.hard_delay_beta,[],condcolors(2,:),'Filled')
% xticklabels({'Easy','Hard'}); title('Delay period pupil beta')
% clean_fig();
% [h,p] = ttest(data.easy_delay_beta,data.hard_delay_beta);
% if p < 0.05
%     plot(1.5,nanmean([nanmean(data.easy_delay_beta), nanmean(data.hard_delay_beta)])+0.2,'*k')
% end
% 
% subplot(2,1,2)
% bar([nanmean(data.easy_stim_beta), nanmean(data.hard_stim_beta)],'w')
% hold on
% scatter(ones(n,1),data.easy_stim_beta,[],condcolors(1,:),'Filled')
% scatter(2.*ones(n,1),data.hard_stim_beta,[],condcolors(2,:),'Filled')
% xticklabels({'Easy','Hard'}); title('Stim presentation beta')
% [h,p] = ttest(data.easy_stim_beta,data.hard_stim_beta);
% if p < 0.05
%     plot(1.5,nanmean([nanmean(data.easy_stim_beta), nanmean(data.hard_stim_beta)])+0.2,'*k')
% end
% clean_fig();
% 
% % Across the board, group-level analyses do not work out
% % Maybe individual analyses will work in TAFKAP precision.
