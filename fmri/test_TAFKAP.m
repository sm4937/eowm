%% Test TAFKAP 
% Run a generate/recover on TAFKAP model fitting 
clear all
addpath(genpath('../../TAFKAP'))

%% Simulate data

nvoxels = 100;
ntrials = 120;
% 180 degrees of stimuli by 75 trials or whatever
stimuli = randperm(180,ntrials);
stim_width = 5;
screen = zeros(ntrials,180);
% Set up noisy stimulus presentation
for trial = 1:ntrials
    screen(trial,:) = normpdf(1:180,stimuli(trial),stim_width);
end
runs = reshape(repmat([1:10],12,1),ntrials,1);

RFs = zeros(50,180);
RF1 = normpdf(1:180,60,15);
RF2 = normpdf(1:180,120,15);
RFs = [repmat(RF1,50,1); repmat(RF2,50,1)];

BOLD = RFs*screen';
% One voxel per row, 1 trial per column
BOLD = BOLD';
% Now one trial per row, 1 voxel per column

test = false(ntrials,1);
test(runs==10) = true;
train = ~test;

p = struct;
p.stimval = stimuli;
p.runNs = runs;
p.test_trials = test;
p.train_trials = train;
p.stim_type = 'circular';

[est, unc, liks, hypers]  = TAFKAP_Decode(BOLD,p)


%% Write some code to actually run TAFKAP, one ROI at a time,
% doing leave-one-run-out cross validation

subject = 4;
% load in BOLD data (already meaned)
load(['PSC/late_delay_mean_PSC_subj' num2str(subject) '_VEselected.mat'])
% laod in task data
load(['task_subj' num2str(subject) '.mat'])
eval(['TAFKAP_subj' num2str(subject) ' = struct; data = task_subj' num2str(subject) ';'])

% If things have been stopped/started, load progress
if contains(ls,['TAFKAP_subj' num2str(subject) '.mat'])
    load(['TAFKAP_subj' num2str(subject) '.mat']);
end

p = struct;
p.runNs = data.overallrun;
p.stim_type = 'circular';
ROI_list = fieldnames(delay_mean_PSC);
ROI_list(contains(ROI_list,{'V2d','V2v','V3d','V3v','VO1','VO2','hV4'})) = [];

for rii = 9:length(ROI_list)
    ROI_name = ROI_list(rii);
    eval(['BOLD = delay_mean_PSC.' ROI_name{1} ';'])
    % load the relevant BOLD data, ROI-specific, where
    % trial is row, and voxel is column
    
    % Make stimulus values go from 0 to 180
    temp = (data.stimval)./2;
    % Trim trials with no stimulus out
    toexclude = isnan(temp) & (sum(isnan(BOLD),2)>0);
    BOLD = BOLD(~toexclude,:); 
    p.stimval = temp(~toexclude);
    task = data(~toexclude,:);
    
    % Trim voxels out with no BOLD (just NaNs)
    BOLD(:,sum(isnan(BOLD),1)>0) = [];
    
    %Z-score each voxel's time course, within each run
    for runii = 1:max(data.overallrun)
        idx = task.overallrun==runii;
        BOLD(idx,:) = zscore(BOLD(idx,:),0,2);
    end
    
    % Cross many-fold fitting, leave one run out
    for fold = 1:max(data.overallrun)
        p.test_trials = task.overallrun == fold;
        p.train_trials = task.overallrun ~= fold;
        % select the 1 run to leave out (to test)

        trial_numbers = find(p.test_trials);
        % determine which rows are being tested, for later saving
        [est, unc, liks, hypers]  = TAFKAP_Decode(BOLD,p);
        
        eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.est(trial_numbers,:) = est;'])
        eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.unc(trial_numbers,:) = unc;'])
        % save results of decoding in subject-specific file
        eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.liks(trial_numbers,:) = liks;'])
        eval(['TAFKAP_subj' num2str(subject) '.' ROI_name{1} '.hypers(trial_numbers,:) = repmat(hypers,length(trial_numbers),1);'])
        % not sure how important this will be later on, but let's save it
        % all
        
        save(['TAFKAP_subj' num2str(subject) '_late_delay.mat'],['TAFKAP_subj' num2str(subject)])
        
    end
end

