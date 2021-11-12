%% Data analysis script for EOWM
% Sarah Master
% 01/2021
clear all
% first-pass at analyzing effort of WM task data
files = dir('data/'); files = char(files.name);
data_folders = files(contains(string(files),'subj'),:);

addpath(genpath('../superbar_all_files'))

% for f = 1:size(data_folders,1) %from existing subj folders, load subject numbers up
%     temp = strsplit(string(data_folders(f,:)),'subj');
%     subjs(f) = str2num(char(temp(2))); 
% end
% subjs(subjs==99) = []; %ignore debug subject
subjs = 4:6; %replace here because pilot 1 is a mess anyway

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
disp(['Between ' num2str(100*min(data.pct_excluded)) '% and ' num2str(100*max(data.pct_excluded)) '% trials excluded per participant.'])
disp(['Mean(std): ' num2str(mean(data.pct_excluded)) '(' num2str(std(data.pct_excluded)) ')']);

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
ax = gca; ax.FontSize = 14;
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
ax = gca; ax.FontSize = 14;

% subplot(2,1,2)
% fulltc = [nanmean([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanmean([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size])];
% full_SEM = [nanstd([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanstd([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size])]./sqrt(n);
% xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)].*(2/1000);
% plot(xs,fulltc(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
% %errorbar(fulltc(2,:),full_SEM(2,:),'Color',condcolors(2,:),'LineWidth',1.5)
% hold on
% plot(xs,fulltc(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
% %errorbar(fulltc(1,:),full_SEM(1,:),'Color',condcolors(1,:),'LineWidth',1.5)
% plot([0 0],ylim,'k--','LineWidth',1.5)
% title('Timecourse of pupil size (correct trials)')
% legend({'Hard','Easy'})
% xlabel('Time (sec)[0 = delay onset]')
% ylabel('Mean pupil size')
% xlim([-0.5 12.5]); ylim([-200 200])
% ax = gca; ax.FontSize = 14;
% fig = gcf; fig.Color = 'w';

subplot(2,2,3)
[h,pval] = ttest(data.cond_early_pupil_size(:,1),data.cond_early_pupil_size(:,2));
superbar([mean(data.cond_early_pupil_size(:,1)) mean(data.cond_early_pupil_size(:,2))],'E',[nanstd(data.cond_early_pupil_size(:,1)) nanstd(data.cond_early_pupil_size(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
ylabel('Mean pupil size 0-6 seconds of delay')
xticklabels(condlabels)
xticks([1:2]);xlim([0.5 2.5])
ax = gca; ax.FontSize = 14;

subplot(2,2,4)
[h,pval] = ttest(data.cond_late_pupil_size(:,1),data.cond_late_pupil_size(:,2));
superbar([mean(data.cond_late_pupil_size(:,1)) mean(data.cond_late_pupil_size(:,2))],'E',[nanstd(data.cond_late_pupil_size(:,1)) nanstd(data.cond_late_pupil_size(:,2))]./sqrt(n),'P',flip(pval*eye(2,2)),'BarFaceColor',condcolors)
ylabel('Mean pupil size 6-12 seconds of delay')
xticklabels(condlabels)
xticks([1:2]);xlim([0.5 2.5])
ax = gca; ax.FontSize = 14;
fig = gcf; fig.Color = 'w';

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
ax = gca; ax.FontSize = 14;

subplot(3,1,2)
fulltc = [nanmean([data.easy_pres_tc_pupil_size data.easy_delay_tc_pupil_size],1);nanmean([data.hard_pres_tc_pupil_size data.hard_delay_tc_pupil_size])];
xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)].*(2/1000);
plot(xs,fulltc(2,:),'Color',condcolors(2,:))
hold on
plot(xs,fulltc(1,:),'Color',condcolors(1,:))
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (correct trials)')
legend({'Hard','Easy'})
xlabel('Time (sec)[0 = delay onset]')
ylabel('Mean pupil size')
xlim([-0.5 12.5]); ylim([-200 200])
ax = gca; ax.FontSize = 14;

subplot(3,1,3)
fulltc = [nanmean([data.easy_pres_tc_pupil_size_incorrect data.easy_delay_tc_pupil_size_incorrect],1);nanmean([data.hard_pres_tc_pupil_size_incorrect data.hard_delay_tc_pupil_size_incorrect])];
xs = [(-length(data.easy_pres_tc_pupil_size_incorrect)+1):0 1:length(data.easy_delay_tc_pupil_size_incorrect)].*(2/1000);
plot(xs,fulltc(2,:),'Color',condcolors(2,:))
hold on
plot(xs,fulltc(1,:),'Color',condcolors(1,:))
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (incorrect trials)')
legend({'Hard','Easy'})
xlabel('Time (sec)[0 = delay onset]')
ylabel('Mean pupil size')
ax = gca; ax.FontSize = 14;



