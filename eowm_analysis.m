%% Data analysis script for EOWM
% Sarah Master
% 01/2021
clear all
% first-pass at analyzing effort of WM task data
files = dir('data/'); files = char(files.name);
data_folders = files(contains(string(files),'subj'),:);
% for f = 1:size(data_folders,1) %from existing subj folders, load subject numbers up
%     temp = strsplit(string(data_folders(f,:)),'subj');
%     subjs(f) = str2num(char(temp(2))); 
% end
% subjs(subjs==99) = []; %ignore debug subject
subjs = 2:4; %replace here because pilot 1 is a mess anyway

data = [];
for s = 1:length(subjs)
    subj = subjs(s);
    data = [data; makeEOWMdatatable(subj)];
end
n = length(data.subj); condlabels = {'Easy','Hard'};

%% Make plots
figure
subplot(2,2,1)
bar(mean(data.cond_accuracy))
hold on
errorbar(mean(data.cond_accuracy),nanstd(data.cond_accuracy)./sqrt(n),'ok')
xticklabels(condlabels)
title('Overall accuracy by condition')

subplot(2,2,2)
bar(nanmean(data.cond_rt))
hold on
errorbar(nanmean(data.cond_rt),nanstd(data.cond_rt)./sqrt(n),'ok')
xticklabels(condlabels)
title('Mean RTs by condition')

subplot(2,2,3)
bar(mean(data.accuracy_by_quad))
hold on
errorbar(mean(data.accuracy_by_quad),nanstd(data.accuracy_by_quad)./sqrt(n),'ok')
title('Accuracy by quadrant')

subplot(2,2,4)
bar(nanmean(data.rt_by_quad))
hold on
errorbar(nanmean(data.rt_by_quad),nanstd(data.rt_by_quad)./sqrt(n),'ok')
title('RT by quadrant')
fig = gcf; fig.Color = 'w';


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
xs = [(-length(data.easy_pres_tc_pupil_size)+1):0 1:length(data.easy_delay_tc_pupil_size)];
plot(xs,fulltc(2,:),'r')
hold on
plot(xs,fulltc(1,:),'b')
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (correct trials)')
legend({'Hard','Easy'})
xlabel('Time (msec)[0 = delay onset]')
ylabel('Mean pupil size')
xlim([-500 6500])
ax = gca; ax.FontSize = 14;

subplot(3,1,3)
fulltc = [nanmean([data.easy_pres_tc_pupil_size_incorrect data.easy_delay_tc_pupil_size_incorrect],1);nanmean([data.hard_pres_tc_pupil_size_incorrect data.hard_delay_tc_pupil_size_incorrect])];
xs = [(-length(data.easy_pres_tc_pupil_size_incorrect)+1):0 1:length(data.easy_delay_tc_pupil_size_incorrect)];
plot(xs,fulltc(2,:),'r')
hold on
plot(xs,fulltc(1,:),'b')
plot([0 0],ylim,'k--','LineWidth',1.5)
title('Timecourse of pupil size (incorrect trials)')
legend({'Hard','Easy'})
xlabel('Time (msec)[0 = delay onset]')
ylabel('Mean pupil size')
ax = gca; ax.FontSize = 14;



