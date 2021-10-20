%% Going through fMRI tools tutorial that Yuna sent from the NYU fMRI Lab course
% Adapting to my purposes as a first pass on fMRI data analysis
clear all
tic 
%specify subject ID as string
subject = 'CC';
sess = 2;
condcolors = [153 0 76; 0 76 153]./255;

datapath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/'];
addpath(genpath(datapath))
addpath(genpath('/Users/sarah/Documents/MATLAB/fmriTools'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/preproc_mFiles'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/vistasoft_ts'))
% pull up some pre-built functions
toc 

% make sure data accessible
savepath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' subject];
newdatadir = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' subject '/fmri/'];
mkdir(newdatadir);
subjfiles = ls([savepath '/fmri']);
ROIpath = [datapath 'ROIs'];

%% Load data up from nifti files

if ~contains(subjfiles,['sess' num2str(sess)]) %fmri data has not already been loaded

    %anatpath = [datapath 'sub-' subject '/ses-anat/anat/sub-' subject '_ses-anat_desc-preproc_T1w.nii'];
    %funcpath = [datapath 'sub-' subject '/ses-func/func/'];
    anatpath = [datapath subject 'anat/'];
    funcpath = [datapath 'sess' num2str(sess) '/'];


    funcfiles = dir(funcpath); runs = [];
    % From subject-specific folder, grab functional run numbers for later use
    for f = 1:length(funcfiles)
        name = funcfiles(f).name;
        if contains(name,'func_volreg_normPctDet')
            file = strsplit(name,'func_volreg_normPctDet');
            runs = [runs; str2double(file{2}(1:2))];
        end
    end

    runs = unique(runs(~isnan(runs)));
    for r = 1:length(runs)
        %functional_tag = ['sub-' subject '_ses-func_task-TASK_run-' num2str(r) '_space-T1w_desc-preproc_bold.nii'];
        if r < 10
            %functional_tag = ['surf_volreg_normPctDet0' num2str(r) '.nii.gz'];
            functional_tag = ['surf_volreg_detrend0' num2str(r) '.nii.gz'];  
        else
            functional_tag = ['surf_volreg_detrend' num2str(r) '.nii.gz'];
        end
        % This is an arbitrary functional file, probably not the one you
        % want to use, really
        filename = [funcpath functional_tag];
        [funcdata{r}] = niftiRead(filename);
    end

    %[anat] = niftiread(anatpath);
    save([savepath '/fmri/sess' num2str(sess) '.mat'],'funcdata','-v7.3');
else
    load([savepath '/fmri/sess' num2str(sess) '.mat']);
end

% Cycle through and display slices, with "slice" command in the dimension
% you want to move over

% figure
% for slice = 1:size(anat,3)
%     imagesc(flipud(squeeze(anat(:,:,slice))'));
%     pause(0.1)
% end
% anat is x, y, z coordinates

% You can go over time with functional data by cycling over the 4th
% dimension
% figure
% ntrs = size(funcdata{1},4); %get number of timepoints
% for t = 1:ntrs
%     imagesc(funcdata{1}(:,:,40,t))
%     drawnow
%     pause(0.01)    
% end

% Check out a single voxel time series
% 
% timeseries = squeeze(funcdata{1}(50,80,50,:));
% figure; plot(timeseries);
% xlabel('time'); ylabel('BOLD');
% fig = gcf; fig.Color = 'w';
% % Cool!

%% So, one thing you could do is match task data to TRs
% One TR = 750 msec
% One trial = ~16 seconds
% Match the timings up to get task-based insights
subject = 4;
datapath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subject)];
allbehavfiles = strsplit(ls(datapath),'.mat'); %get behavioral data to match up

TR = 0.75; % in seconds

long = [];
for sess = 1:2
    
    load([savepath '/fmri/sess' num2str(sess) '.mat']);
    % load up functional data for that session
    runs = length(funcdata);
    nTRs = size(funcdata{1}.data,4); %how many total TRs?
    
    for r = 1:runs
        run = r-1; long_run = [];
        behavidx = find(contains(allbehavfiles,['run' num2str(run) '_sess' num2str(sess)]));
        % find and load the behav file for this run
        beginidx = find(allbehavfiles{behavidx}=='o'); beginidx = beginidx(1)-1;
        behav = load([datapath '/' allbehavfiles{behavidx}(beginidx:end) '.mat']);
        times = 0:TR:(TR*nTRs); times = times(1:nTRs)';
        run_start = behav.p.trial_start(1);
        run_end = behav.p.iti_start(end)-run_start;

        all_times = [behav.p.trial_start behav.p.targ_start behav.p.delay_start behav.p.test_start behav.p.feedback_start behav.p.iti_start behav.p.iti_start+behav.p.itis];
        % make a column vector which starts at 0 and ends at the end of trial
        all_times = reshape([all_times-run_start]',numel(all_times),1);
        all_times = [0; behav.p.start_wait+all_times; all_times(end)+behav.p.end_wait];

        
        % mark all events that occur
        % 1 is trial start, 2 is target start, 3 is delay start
        % 4 is test start, 5 is feedback start, 6 is iti start, 7 is 
        % time betwen ITI end and next trial start
        all_events = [0; repmat([1; 2; 3; 4; 5; 6; 7],12,1); 0];
        all_conds = [0; reshape(repmat(behav.p.conditions(:,1),1,7)',84,1); 0];
        all_ang = [NaN; reshape(repmat(behav.p.wm_ang,1,7)',84,1); NaN];
        all_trials = repmat([1:12],7,1); all_trials = [0; all_trials(:); 0];
        % 0 padding and extra seconds added on come from wait time at start
        % and end of each run
        
        for tt = 1:length(times)
            slice = times(tt);
            idx = find(all_times<slice);
            if ~isempty(idx)
                long_run = [long_run; subject sess r all_trials(idx(end)) slice all_events(idx(end)) all_conds(idx(end)) all_ang(idx(end))];
            end
        end
        
        % get relative time in delay period to compare late delay to early
        % delay later on, etc.
        delay_times = zeros(size(long_run,1),1);
        for trial = 1:max(unique(long_run(:,4)))
            start = behav.p.start_wait + behav.p.delay_start(trial) - run_start;
            delay_times(long_run(:,6)==3&long_run(:,4)==trial,:) = long_run(long_run(:,6)==3&long_run(:,4)==trial,5) - start
        end
        long_run = [long_run delay_times];
        long = [long; long_run];
        
    end
end
% label in table for clarity
full = table; 
full.subj = long(:,1);
full.sess = long(:,2);
full.run = long(:,3); 
full.trial = long(:,4);
full.time = long(:,5);
full.epoch = long(:,6); 
full.cond = long(:,7);
full.stimval = long(:,8);
full.delay_times = long(:,9);
full.overalltrial = cumsum([0; diff(full.trial)~=0]);
full.overallrun = cumsum([1; diff(full.run)~=0]);

% delay_times = zeros(height(full),1);
% full.delay_times = delay_times;

idxes = ([false; diff(full.overalltrial)~=0]);
eval(['task_subj' num2str(subject) ' = full(idxes,:);'])
save(['task_subj' num2str(subject) '.mat'], ['task_subj' num2str(subject)])

%% Create GLM 
% Predictors are each epoch, each condition

tau = 2; %tau/10 is rate of decay
delta = 5; %experiment-derived time lag parameter

% Set the HRF according to these parameters
t = [0:1:30];
tshift = max(t-delta,0);
HIRF = (tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);

predictors(:,1) = ones(height(full),1); %baseline image intensity (intercept term)
predictors(:,2) = repmat([1:(nTRs-1)]',runs*sess,1); %linear drift component
% in full.event:
% 1 is trial start, 2 is target start, 3 is delay start
% 4 is test start, 5 is feedback start, 6 is iti start, 7 is 
% time betwen ITI end and next trial start
predictors(:,3) = full.epoch==2 & full.cond==1; %target stimulus on screen
predictors(:,4) = full.epoch==3 & full.cond==1; %delay period
predictors(:,5) = full.epoch==6 & full.cond==1; %ITI
% redo all the above predictors for hard trials
predictors(:,6) = full.epoch==2 & full.cond==2; %target stimulus on screen
predictors(:,7) = full.epoch==3 & full.cond==2; %delay period
predictors(:,8) = full.epoch==6 & full.cond==2; %ITI

dummies = predictors;
dummies(:,3:8) = dummies(:,3:8) .* [3:8];

model = pinv(predictors);

voxeltc = funcdata{r}(x,y,z,:); 
%grab timecourse (tc) of each voxel
betas = model * squeeze(voxeltc);


%% Analyze functional data in terms of percent signal change,
% across all PRF-defined ROIs 

subject = 'CC';
subjnum = 4;

ROI_folder = dir(ROIpath);
ROI_list = char(ROI_folder.name); %list all possible ROIs
ROI_list(1:2,:) = []; % get rid of blank indices in folder
%ROI_list = string(ROI_list);

%PRFparams = niftiread('/System/Volumes/Data/d/DATA/data/vRF_tcs/CC/RF1/CC_RF1_vista/RF_surf_25mm-fFit.nii.gz');
%save(['PRFparams_subj' num2str(subjnum) '.mat'],'PRFparams')
load(['PRFparams_subj' num2str(subjnum) '.mat'])
% Params is 8D, and each dimension is:
% pol : polar angle
% ve : variance explained
% ecc : eccentricity
% sigmamajor : sigma of the gaussian describing the RF
% exponent : the exponent of the gaussian RF
% x0 : center of the PRF in x coords
% y0 : center of the PRF in y coords
% b : not sure

hemi = 'bilat';
% pick the hemisphere you want to examine
hemi_list = ROI_list(contains(string(ROI_list),hemi),:);

%Pre-define some variables that you want averages over
delay_mean_PSC = struct();

for area = 1:size(hemi_list,1)
    ROI_filename = strrep(hemi_list(area,:),' ','');
    ROI = niftiRead([ROIpath '/' ROI_filename]);
    pat = '.(\w*).nii.gz';
    ROI_name = regexp(ROI_filename, pat, 'tokens');
    trial_index = [];
    all_trial_signal = [];
    
    % get variance explained from PRF fitting
    VE = PRFparams(:,:,:,2); VE = VE(ROI.data>0);
    nvoxels = sum(VE>0.1);
    
    %initialize saving structure here
    eval(['delay_mean_PSC.' + string(ROI_name) + ' = NaN(max(full.overalltrial),nvoxels);'])

    for sess = 1:2
        load([savepath '/fmri/sess' num2str(sess) '.mat']);
        % load up functional data for that session
        runs = length(funcdata);

        for r = 1:runs %cycle over runs within that session
            
            all_voxel_tcs = niftiExtract(funcdata{r},ROI);
            % Tommy wrote this function to extract ROI-based functional
            % maps
            
            all_voxel_tcs = all_voxel_tcs';
            spec_voxel_tcs = all_voxel_tcs(VE>0.1,:);
            
            % Z-score of this?
                         
            % reshape the voxels into a 2D matrix,
            % where time is the 2nd dimension and voxel # is the first 
            nTRs = size(funcdata{r}.data,4);

            fmriSignal = spec_voxel_tcs;
            fmriResponse = 100 * ((fmriSignal./(nanmean(fmriSignal,2))) - 1);
            %fmriResponse = fmriSignal;
            
            data = full(full.run==r&full.sess==sess,:);
            % Cycle over trials 
            for t = 1:max(unique(data.trial))
                trial_signal = NaN(nvoxels,55);
                trial_signal(:,1:sum(data.trial==t)) = fmriResponse(:,data.trial==t);
                % Grab TR's inside each trial, in order

                all_trial_signal = [all_trial_signal; trial_signal];
                % Grab percent signal change in that voxel over each of those
                % trials
                overalltrial = unique(data.overalltrial(data.trial==t));
                trial_index = [trial_index; repmat(overalltrial,size(trial_signal,1),1)];
                % Keep an index of which trial is in which row (remember that
                % each row is a voxel's time course inside that trial, and that
                % we're averaging across voxels, over all trials. In the
                % future you may want to use this index to compare percent
                % signal change across trial types.)
                
                delay = (data.epoch(data.trial==t)==3)';
                late_delay = data.delay_times(data.trial==t) >= 6;
                eval(['delay_mean_PSC.' + string(ROI_name) + "(overalltrial,:) = nanmean(trial_signal(:,late_delay),2)';"])
                % save mean PSC in delay period on this trial, for that ROI
                % resulting struct has trial (rows) x voxel (columns)
                % measurements of mean delay period PSC
            end % of cycling over trials

        end % of cycling over runs
    
    end % of cycling over sessions
    
    figure(5)
    hard_trials = unique(full.overalltrial(full.cond==2,:));
    hard_idx = sum(trial_index == hard_trials',2)==1;
    errorbar(nanmean(all_trial_signal(hard_idx,:)),nanstd(all_trial_signal(hard_idx,:))./sqrt(nvoxels*sess),'Color',condcolors(1,:),'LineWidth',1.25,'DisplayName','Hard trials')
    xticks(1:22)
    xticklabels(0:TR:TR*22)
    hold on
    errorbar(nanmean(all_trial_signal(~hard_idx,:)),nanstd(all_trial_signal(~hard_idx,:))./sqrt(nvoxels*sess),'Color',condcolors(2,:),'LineWidth',1.25,'DisplayName','Easy trials')
    plot([2.5 2.5],ylim,'k--','LineWidth',1.25,'DisplayName','Stim onset')
    plot([20 20],ylim,'k--','LineWidth',1.25,'DisplayName','Test onset')
    plot([21.5 21.5],ylim,'k--','LineWidth',1.25,'DisplayName','Feedback onset')
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 12;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean percent signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    
    pause(5); saveas(fig,['PSC/' + string(ROI_name{1}) + '.jpg']);
    close 5
    
    save(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected.mat'],'delay_mean_PSC')
end % of cycling over areas of interest (ROIs)

%% %% Convolving with an HRF
% each presentation of something on screen, like the target & test stimuli
% and the feedback, could be thought of as producing a response
% but the measurement of that response is complicated by how slow the HRF
% actually is
% approximate how the neural activity from these events turns into fMRI by
% convolving these events with an approximate HRF
r = 1;
neuralActivity = full.event(full.run==r);

% time constant, time offset parameters
tau = 2; %tau/10 is rate of decay
delta = 5; %experiment-derived time lag parameter

% Plot the HIRF with these parameter values:
t = [0:1:30];
tshift = max(t-delta,0);
HIRF = (tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);
% figure(1); clf;
% % Plot it
% plot(t,HIRF);
% title('Hemodynamic Impulse Response Function')
% ylabel('Hemodynamic response')
% xlabel('Time (sec)')

% Now apply the HIRF to the neuralActivity and see what happens
fmriSignal = conv(neuralActivity,HIRF);
noise = 0.1 * randn(size(fmriSignal));
drift = 0.01 * [1:length(fmriSignal)];

fmriSignal = fmriSignal + noise + drift; 
figure
plot(fmriSignal(1:200))
hold on
plot(neuralActivity(1:200))
xlabel('Time'); ylabel('Signal intensity')

fmriResponse = 100 * ((fmriSignal(1:length(neuralActivity))...
    /(mean(fmriSignal(1:length(neuralActivity)))) - 1));
%figure; plot(fmriResponse)

% this calculates percent signal change compared to baseline, which is the
% mean of the overall signal (without the 0 padding that Matlab's conv
% tacks on to the end)

%% Running a GLM
% Create the design matrix, X, for the GLM

%create model components
predictors(:,1) = fmriSignal(1:length(neuralActivity))'; % HRF component
predictors(:,2) = [1:length(neuralActivity)]'; %linear drift component
predictors(:,3) = ones(length(neuralActivity),1); %baseline image intensity (intercept term)

X = size(funcdata{r},1); Y = size(funcdata{r},2);
Z = size(funcdata{r},3); nTRs = length(neuralActivity);

model = pinv(predictors); betas = NaN(X,Y,Z,size(model,1));
for x = 1:X %cycle over X coordinates in brain space
    for y = 1:Y
        for z = 1:Z
            voxeltc = funcdata{r}(x,y,z,:); 
            %grab timecourse (tc) of each voxel
            betas(x,y,z,:) = model * squeeze(voxeltc);
        end
    end
end

% Plot betas for each voxel
figure
subplot(3,1,1)
betas_neural = betas(:,:,:,1);
plot(betas_neural(:))
title('Estimates of neural activity')
ylabel('Amplitude (AU)')
xlabel('Voxel')

subplot(3,1,2)
betas_drift = betas(:,:,:,2);
plot(betas_drift(:))
title('Estimates of drift')
ylabel('Amplitude (AU)')
xlabel('Voxel')

subplot(3,1,3)
betas_baseline = betas(:,:,:,3);
plot(betas_baseline(:))
title('Estimates of baseline')
ylabel('Amplitude (AU)')
xlabel('Voxel')
fig = gcf; fig.Color = 'w';


