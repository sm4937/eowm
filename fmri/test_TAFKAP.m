%% Test TAFKAP 
% Run a generate/recover on TAFKAP model fitting 
clear all
addpath(genpath('../../TAFKAP'))

% Simulate a fake experiment, fake data, fake RFs
% Decode these stimuli to see how it does (it's perfect)

nvoxels = 350;
ntrials = 240;
% 180 degrees of stimuli by 75 trials or whatever
stimuli = randsample(360,ntrials,true);
stim_width = 5;

% Set up noisy stimulus presentation
screen = normpdf(1:360,stimuli,stim_width);

runs = reshape(repmat([1:20],12,1),ntrials,1);

% RF1 = normpdf(1:180,60,15);
% RF2 = normpdf(1:180,120,15);
% RFs = [repmat(RF1,50,1); repmat(RF2,50,1)];

RF_stds = normrnd(6,10,nvoxels,1);
RF_stds(RF_stds<1) = 2;
RF_means = randsample(360,nvoxels,true);
RFs = normpdf(1:360,RF_means,RF_stds);

BOLD = RFs*screen';
% One voxel per row, 1 trial per column
BOLD = BOLD';
% Now one trial per row, 1 voxel per column

for fold = 1:max(runs)

    test = false(ntrials,1);
    test(runs==fold) = true;
    train = ~test;

    p = struct;
    p.stimval = stimuli./2;
    p.runNs = runs;
    p.test_trials = test;
    p.train_trials = train;
    p.stim_type = 'circular';
    p.Nboot = 1000;

    [est(test,:), unc(test,:), liks(test,:), hypers(fold,:)]  = TAFKAP_Decode(BOLD,p);
end

figure; subplot(1,2,1)
scatter(stimuli,est,'Filled')
title('Stimulus estimates versus real stimuli: Generate/Recover')
ylabel('Estimates')
xlabel('Real')

subplot(1,2,2)
histogram(unc,'DisplayName','Estimation')
title('Estimated uncertainty')
xlabel('Estimated unc')
ylabel('Frequency')
hold on
plot([6 6], ylim, '--k','DisplayName','Mean uncertainty of channels')
fig = gcf; fig.Color ='w';
legend('location','best')

%% Now write another gen/rec, but for this one, use your real voxel properties
% and N's from your experiment/task
clear all

subjnum = 4;
subject = 'CC';

load(['PRFparams_subj' num2str(subjnum) '.mat'])
load(['task_subj' num2str(subjnum) '.mat'])
eval(['task = task_subj' num2str(subjnum) ';'])

ROI_name = 'IPS0';
% supposed to be a very precise memory area, despite being small in nvoxels
ROI = niftiRead(['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/' ...
 'ROIs/bilat.' ROI_name '.nii.gz']);

% get variance explained from PRF fitting, real nvoxels to use
VE = PRFparams(:,:,:,2); VE = VE(ROI.data>0);
nvoxels = sum(VE>0.1);

voxel_RF_centers = rad2deg(PRFparams(:,:,:,1)); 
voxel_RF_centers = voxel_RF_centers(ROI.data>0);
voxel_RF_sigmas = PRFparams(:,:,:,4);
voxel_RF_sigmas = voxel_RF_sigmas(ROI.data>0);

RFs = normpdf(0:359,voxel_RF_centers,voxel_RF_sigmas);
RFs = RFs(VE>0.1,:);

% Get real ntrials, real stimulus values
ntrials = sum(~isnan(task.stimval));
runs = task.overallrun(~isnan(task.stimval));
stimuli = task.stimval(~isnan(task.stimval));

stim_width = 3;

% Set up noisy stimulus presentation
screen = normpdf(0:359,stimuli,stim_width);

BOLD = RFs*screen';
% One voxel per row, 1 trial per column
BOLD = BOLD';
% Now one trial per row, 1 voxel per column
BOLD = BOLD + normrnd(0,1,size(BOLD)) + normrnd(0,1,size(BOLD,1),1); %add a tiny amount of normally distributed noise
%BOLD = zscore(BOLD,[],1);

for fold = 1:max(runs)

    test = false(ntrials,1);
    %test(sum(runs==[fold:fold+4],2)==1) = true;
    test(runs==fold) = true;
    train = ~test;

    p = struct;
    p.stimval = stimuli./2;
    p.runNs = runs;
    p.test_trials = test;
    p.train_trials = train;
    p.stim_type = 'circular';
    p.Nboot = 1000;

    [est(test,:), unc(test,:), liks(test,:), hypers(fold,:)]  = TAFKAP_Decode(BOLD,p);
end

figure; subplot(1,2,1)
scatter(stimuli,est.*2,'Filled')
title('Stimulus estimates versus real stimuli: Generate/Recover')
ylabel('Estimates')
xlabel('Real')

subplot(1,2,2)
histogram(unc,'DisplayName','Estimation')
title('Estimated uncertainty')
xlabel('Estimated unc')
ylabel('Frequency')
hold on
plot([6 6], ylim, '--k','DisplayName','Mean uncertainty of channels')
fig = gcf; fig.Color ='w';
legend('location','best')



%% Write some code to actually run TAFKAP, one ROI at a time,
% doing leave-one-run-out cross validation
% clear all
% 
% subject = 6;
% % load in BOLD data (already meaned)
% load(['PSC/late_delay_mean_PSC_subj' num2str(subject) '_VEselected.mat'])
% % laod in task data
% load(['task_subj' num2str(subject) '.mat'])
% eval(['TAFKAP_subj' num2str(subject) ' = struct; data = task_subj' num2str(subject) ';'])
% 
% % If things have been stopped/started, load progress
% if contains(ls,['TAFKAP_subj' num2str(subject) '_late_delay.mat'])
%     load(['TAFKAP_subj' num2str(subject) '_late_delay.mat']);
% end
% 
% p = struct;
% p.runNs = data.overallrun;
% p.stim_type = 'circular';
% %p.Nboot = 1000;
% ROI_list = fieldnames(delay_mean_PSC);
% ROI_list(contains(ROI_list,{'V2d','V2v','V3d','V3v','VO1','VO2','LO1','LO2','TO1','TO2'})) = [];
% 
% figure;
% for rii = 1:length(ROI_list)
%     ROI_name = ROI_list(rii);
%     eval(['BOLD = delay_mean_PSC.' ROI_name{1} ';'])
%     % load the relevant BOLD data, ROI-specific, where
%     % trial is row, and voxel is column
%     
%     % Make stimulus values go from 0 to 180
%     temp = (data.stimval)./2;
%     % Trim trials with no stimulus out
%     toexclude = isnan(temp) & (sum(isnan(BOLD),2)>0);
%     BOLD = BOLD(~toexclude,:); 
%     p.stimval = temp(~toexclude);
%     task = data(~toexclude,:);
%     
%     % Trim voxels out with no BOLD (just NaNs)
%     BOLD(:,sum(isnan(BOLD),1)>0) = [];
%     
%     %Z-score each voxel's time course, within each run
%     for runii = 1:max(data.overallrun)
%         idx = task.overallrun==runii;
%         BOLD(idx,:) = zscore(BOLD(idx,:),0,1);
%     end
%     
%     % Cross many-fold fitting, leave one run out
%     for fold = 1:max(data.overallrun)
%         p.test_trials = task.overallrun == fold;
%         p.train_trials = task.overallrun ~= fold;
%         % select the 1 run to leave out (to test)
% 
%         trial_numbers = find(p.test_trials);
%         % determine which rows are being tested, for later saving
%         [est, unc, liks, hypers]  = TAFKAP_Decode(BOLD,p);
%         
%         eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.est(trial_numbers,:) = est;'])
%         eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.unc(trial_numbers,:) = unc;'])
%         % save results of decoding in subject-specific file
%         eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.liks(trial_numbers,:) = liks;'])
%         eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.hypers(trial_numbers,:) = repmat(hypers,length(trial_numbers),1);'])
%         % not sure how important this will be later on, but let's save it
%         % all
%         
%         save(['TAFKAP_subj' num2str(subject) '_late_delay.mat'],['TAFKAP_subj' num2str(subject)])
%         
%     end
%     
%     eval(['estimates = TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.est;'])
%     
%     subplot(3,4,rii)
%     scatter(p.stimval(~isnan(p.stimval)),estimates(~isnan(p.stimval)),'Filled')
%     fig = gcf; fig.Color = 'w';
%     title(['Decoding accuracy of ' ROI_name{1}])
%     
% end

subject = 4;
% Grab task data
load(['task_subj' num2str(subject) '.mat'])
eval(['task = task_subj' num2str(subject) ';'])

%Specify which file you're interested in
spec = '_late_delay';

if contains(ls,['TAFKAP_subj' num2str(subject) spec '.mat'])
    % if TAFKAP has been run on this subject already, load that output
    load(['TAFKAP_subj' num2str(subject) spec '.mat'], 'TAFKAP_output')
    
else %TAFKAP has not been run on this subject yet
    % load in BOLD data (already meaned)
    load(['PSC/late_delay_mean_PSC_subj' num2str(subject) '_VEselected.mat'])

    BOLD = delay_mean_PSC;
    TAFKAP_output = run_TAFKAP_eowm(BOLD,task);
    save(['TAFKAP_subj' num2str(subject) spec '.mat'], 'TAFKAP_output')
    
end % of deciding to run TAFKAP or not


ROI_list = fieldnames(TAFKAP_output);

figure
% Cycle through each ROI and plot estimation accuracy
for rii = 1:length(ROI_list)
    ROI_name = ROI_list{rii};
    
    eval(['estimates = TAFKAP_output.' ROI_name '.est;'])
    
    subplot(3,4,rii)
    scatter(task.stimval(~isnan(task.stimval)),estimates,'Filled')
    fig = gcf; fig.Color = 'w';
    title(['Decoding accuracy of ' ROI_name])
    
end

