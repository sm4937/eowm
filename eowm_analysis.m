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

disp(['Mean delay period pupil size easy: ' num2str(mean(data.cleaned_pupilsize_bycondition(:,1)))])
disp(['Mean delay period pupil size hard: ' num2str(mean(data.cleaned_pupilsize_bycondition(:,2)))])

disp(['Between ' num2str(100*min(data.pct_excluded)) '% and ' num2str(100*max(data.pct_excluded)) '% trials excluded per participant.'])
disp(['Mean(std): ' num2str(mean(data.pct_excluded)) '(' num2str(std(data.pct_excluded)) ')']);

pupil.subj = subjs';
pupil.pres = [data.alltrials_pres_mean_pupil_size];
pupil.earlydelay = data.alltrials_earlydelay_mean_pupil_size;
pupil.delay = data.alltrials_delay_mean_pupil_size;
pupil.latedelay = data.alltrials_latedelay_mean_pupil_size;
save('data/pupil.mat','pupil')

%% Make plots
% Condition accuracy, with t-tests
figure
subplot(2,2,1)
% errorbar(1, mean(data.cond_accuracy(:,1)),nanstd(data.cond_accuracy(:,1))./sqrt(n),'o','LineWidth',2,'Color',condcolors(1,:),'DisplayName','Easy trials')
% hold on
% errorbar(2, mean(data.cond_accuracy(:,2)),nanstd(data.cond_accuracy(:,2))./sqrt(n),'o','LineWidth',2,'Color',condcolors(2,:),'DisplayName','Hard trials')
[h,pval] = ttest(data.cond_accuracy(:,1),data.cond_accuracy(:,2));
superbar([mean(data.cond_accuracy(:,1)) mean(data.cond_accuracy(:,2))],'E',[nanstd(data.cond_accuracy(:,1)) nanstd(data.cond_accuracy(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
xticklabels(condlabels)
xticks([1:2]); xlim([0.5 2.5]); 
ax = gca; ax.FontSize = 18;
ylabel('% Correct');

subplot(2,2,2)
% errorbar(1, mean(data.cond_rt(:,1)),nanstd(data.cond_rt(:,1))./sqrt(n),'o','LineWidth',2,'Color',condcolors(1,:),'DisplayName','Easy trials')
% hold on
% errorbar(2, mean(data.cond_rt(:,2)),nanstd(data.cond_rt(:,2))./sqrt(n),'o','LineWidth',2,'Color',condcolors(2,:),'DisplayName','Hard trials')
[h,pval] = ttest(data.cond_rt(:,1),data.cond_rt(:,2));
superbar([mean(data.cond_rt(:,1)) mean(data.cond_rt(:,2))],'E',[nanstd(data.cond_rt(:,1)) nanstd(data.cond_rt(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
xticklabels(condlabels)
xticks([1:2]);xlim([0.5 2.5])
ylabel('Mean RT (seconds)')
ax = gca; ax.FontSize = 18;

subplot(2,1,2)

fulltc = [nanmean([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanmean([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)];
full_SEM = [nanstd([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanstd([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)]./sqrt(n);
toclean = sum(isnan(fulltc),1)>0 & sum(isnan(full_SEM),1)>0;
fulltc(:,toclean) = []; full_SEM(:,toclean) = [];

ntrim = 10;
% trim last n points from full TC, they're very messy
fulltc = fulltc(:,1:end-ntrim);
full_SEM = full_SEM(:,1:end-ntrim);

xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)].*(2/1000);
xs(:,toclean) = [];
xs = xs(:,1:end-ntrim);
%errorbar(xs,fulltc(2,:),full_SEM(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
%plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
%errorbar(xs,fulltc(1,:),full_SEM(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
%plot([0 0],ylim,'k--','LineWidth',1.5)
btwn_fill = [fulltc(2,:) + full_SEM(2,:), fliplr(fulltc(2,:))-fliplr(full_SEM(2,:))];
fill_xs = [xs, fliplr(xs)];
fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials');
hold on
btwn_fill = [fulltc(1,:) + full_SEM(1,:), fliplr(fulltc(1,:))-fliplr(full_SEM(1,:))];
fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials');
plot([0 0],ylim,'k--','LineWidth',1.5)
plot(xs,fulltc(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
title('Timecourse of pupil size (correct trials)')
legend({'Hard','Easy'})
xlabel('Time (seconds) [0 = delay onset]')
ylabel('Mean % signal change: pupil size')
[fig,ax] = clean_fig();

% subplot(2,2,3)
% [h,pval] = ttest(data.cond_early_pupil_size(:,1),data.cond_early_pupil_size(:,2));
% superbar([mean(data.cond_early_pupil_size(:,1)) mean(data.cond_early_pupil_size(:,2))],'E',[nanstd(data.cond_early_pupil_size(:,1)) nanstd(data.cond_early_pupil_size(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
% ylabel('Mean pupil size 0-6 seconds of delay')
% xticklabels(condlabels)
% xticks([1:2]);xlim([0.5 2.5])
% ax = gca; ax.FontSize = 14;
% 
% subplot(2,2,4)
% [h,pval] = ttest(data.cond_late_pupil_size(:,1),data.cond_late_pupil_size(:,2));
% superbar([mean(data.cond_late_pupil_size(:,1)) mean(data.cond_late_pupil_size(:,2))],'E',[nanstd(data.cond_late_pupil_size(:,1)) nanstd(data.cond_late_pupil_size(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
% ylabel('Mean pupil size 6-12 seconds of delay')
% xticklabels(condlabels)
% xticks([1:2]);xlim([0.5 2.5])
% ax = gca; ax.FontSize = 14;
% fig = gcf; fig.Color = 'w';

% subplot(2,2,3)
% bar(mean(data.accuracy_by_quad))
% hold on
% errorbar(mean(data.accuracy_by_quad),nanstd(data.accuracy_by_quad)./sqrt(n),'ok')
% title('Accuracy by quadrant')
% 
% subplot(2,2,4)
% bar(nanmean(data.rt_by_quad))
% hold on
% errorbar(nanmean(data.rt_by_quad),nanstd(data.rt_by_quad)./sqrt(n),'ok')
% title('RT by quadrant')
% fig = gcf; fig.Color = 'w';


% display pupil sizes in two conditions
[h,p1] = ttest([data.cond_pres_pupil_size(:,2)-data.cond_pres_pupil_size(:,1)]);
[h,p2] = ttest([data.cond_early_pupil_size(:,2)-data.cond_early_pupil_size(:,1)]);
[h,p3] = ttest([data.cond_delay_pupil_size(:,2)-data.cond_delay_pupil_size(:,1)]);
[h,p4] = ttest([data.cond_late_pupil_size(:,2)-data.cond_late_pupil_size(:,1)]); %is the difference different than 0?

figure
subplot(3,1,1)
bar([nanmean(data.cond_pres_pupil_size(:,2)-data.cond_pres_pupil_size(:,1)) nanmean(data.cond_early_pupil_size(:,2)-data.cond_early_pupil_size(:,1)) nanmean(data.cond_delay_pupil_size(:,2)-data.cond_delay_pupil_size(:,1)) nanmean(data.cond_late_pupil_size(:,2)-data.cond_late_pupil_size(:,1))])
hold on
errorbar([nanmean(data.cond_pres_pupil_size(:,2)-data.cond_pres_pupil_size(:,1)) nanmean(data.cond_early_pupil_size(:,2)-data.cond_early_pupil_size(:,1)) nanmean(data.cond_delay_pupil_size(:,2)-data.cond_delay_pupil_size(:,1)) nanmean(data.cond_late_pupil_size(:,2)-data.cond_late_pupil_size(:,1))], ...
    [nanstd(data.cond_pres_pupil_size(:,2)-data.cond_pres_pupil_size(:,1)) nanstd(data.cond_early_pupil_size(:,2)-data.cond_early_pupil_size(:,1)) nanstd(data.cond_delay_pupil_size(:,2)-data.cond_delay_pupil_size(:,1)) nanstd(data.cond_late_pupil_size(:,2)-data.cond_late_pupil_size(:,1))]./sqrt(n),'ok')
title('Pupil size by epoch')
ylabel('Hard condition pupil size - easy condition size')
xticklabels({'Stim presentation','Early delay','Whole delay period','Late delay'})
fig = gcf; fig.Color = 'w';
for cond = 1:4
    eval(['p = p' num2str(cond) ';'])
    symbol = ['k' get_sig_symbol(p)];
    plot(cond,150,symbol)
end
xtickangle(30)
ax = gca; ax.FontSize = 18;

subplot(3,1,2)
fulltc = [nanmean([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanmean([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)];
full_SEM = [nanstd([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanstd([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)]/sqrt(n*120);
xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)].*(2/1000);
errorbar(xs,fulltc(2,:),full_SEM(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
%plot(xs,fulltc(2,:),'Color',condcolors(2,:))
hold on
%plot(xs,fulltc(1,:),'Color',condcolors(1,:))
errorbar(xs,fulltc(1,:),full_SEM(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (correct trials)')
legend({'Hard','Easy'})
xlabel('Time (sec)[0 = delay onset]')
ylabel('Mean pupil size')
ax = gca; ax.FontSize = 18;

subplot(3,1,3)
fulltc = [nanmean([data.easy_pres_tc_pupil_size_incorrect data.easy_delay_tc_pupil_size_incorrect],1);nanmean([data.hard_pres_tc_pupil_size_incorrect data.hard_delay_tc_pupil_size_incorrect],1)];
full_SEM = [nanstd([data.easy_pres_tc_pupil_size_incorrect data.easy_delay_tc_pupil_size_incorrect],1);nanstd([data.hard_pres_tc_pupil_size_incorrect data.hard_delay_tc_pupil_size_incorrect],1)]/sqrt(n*120);
xs = [(-length(data.easy_pres_tc_pupil_size_incorrect)+1):0 1:length(data.easy_delay_tc_pupil_size_incorrect)].*(2/1000);
errorbar(xs,fulltc(2,:),full_SEM(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
hold on
errorbar(xs,fulltc(1,:),full_SEM(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (incorrect trials)')
legend({'Hard','Easy'})
xlabel('Time (sec)[0 = delay onset]')
ylabel('Mean % signal change: pupil size')
ax = gca; ax.FontSize = 18;

%% Make plots appropriate for VSS Poster 2022

% Condition accuracy, withOUT t-tests
% Not averaged across subjects, but showing all subjects' behavior


figure
subplot(2,2,1)
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
plot(xlim,[0.7 0.7],'-.','LineWidth',1,'Color','k','DisplayName','Target accuracy for hard trials')
plot(xlim,[0.9 0.9],'-.','LineWidth',1,'Color','k','DisplayName','Target accuracy for easy trials')
[fig,ax] = clean_fig();
ylabel('% correct');

subplot(2,2,2)
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

subplot(2,1,2)
fulltc = [nanmean([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanmean([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)];
full_SEM = [nanstd([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanstd([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size],1)]./sqrt(n);
toclean = sum(isnan(fulltc),1)>0 & sum(isnan(full_SEM),1)>0;
fulltc(:,toclean) = []; full_SEM(:,toclean) = [];

ntrim = 10;
% trim last n points from full TC, they're very messy
fulltc = fulltc(:,1:end-ntrim);
full_SEM = full_SEM(:,1:end-ntrim);

xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)].*(2/1000);
xs(:,toclean) = [];
xs = xs(:,1:end-ntrim);
%errorbar(xs,fulltc(2,:),full_SEM(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
%plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
%errorbar(xs,fulltc(1,:),full_SEM(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
%plot([0 0],ylim,'k--','LineWidth',1.5)
btwn_fill = [fulltc(2,:) + full_SEM(2,:), fliplr(fulltc(2,:))-fliplr(full_SEM(2,:))];
fill_xs = [xs, fliplr(xs)];
fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials');
hold on
btwn_fill = [fulltc(1,:) + full_SEM(1,:), fliplr(fulltc(1,:))-fliplr(full_SEM(1,:))];
fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials');
plot(xs,fulltc(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
title('Timecourse of pupil size (correct trials)')
xlim([-0.5 12])
plot([0 0],ylim,'k--','LineWidth',1.5)
legend({'Hard','Easy'})
xlabel('Time (seconds) [0 = delay onset]')
ylabel('Mean % signal change: pupil size')
[fig,ax] = clean_fig();

[h,p] = ttest2(data.cleaned_pupilsize_bycondition(:,1),data.cleaned_pupilsize_bycondition(:,2));
figure; errorbar(1:2,nanmean(data.cleaned_pupilsize_bycondition),nanstd(data.cleaned_pupilsize_bycondition)/sqrt(n),'k','LineWidth',1.5)


%% Let's do a new pupil size analysis - cleaner, more principled
% In the style of Wiehler et al., 2022


