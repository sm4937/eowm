%% Going through fMRI tools tutorial that Yuna sent from the NYU fMRI Lab course
% Adapting to my purposes as a first pass on fMRI data analysis
clear all

addpath(genpath('/Users/sarah/Documents/MATLAB/fmriTools'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/preproc_mFiles'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/vistasoft_ts'))
% pull up some pre-built functions

condcolors = [153 0 76; 0 76 153]./255;

%% Load data up from nifti files

%specify subject ID as string
load('subjinits.mat')

subjnum = 4;
subject = subjinits{subjnum};
sess = 2;

% make sure data accessible
savepath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];
newdatadir = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum) '/fmri/'];
datapath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/'];
mkdir(newdatadir);
subjfiles = ls(newdatadir);
ROIpath = [datapath 'ROIs'];

if ~contains(subjfiles,['sess' num2str(sess)]) %fmri data has not already been loaded
    pull_fmri_data(subjnum,sess)
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
% One trial = ~20 seconds
% Match the timings up to get task-based insights

full = get_full_table(subjnum,sess);

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

% Get task data
load(['task_subj' num2str(subjnum) '.mat'])
eval(['task = task_subj' num2str(subjnum) ';'])

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
% y0 : center of the PRF in y coords, flipped from what you would think
% i.e. negative number are above the horizontal meridian
% b : amplitude

hemi = 'bilat';
% pick the hemisphere you want to examine
hemi_list = ROI_list(contains(string(ROI_list),hemi),:);
% trim ROIs which not all subjects have, for comparison's sake
hemi_list(contains(string(hemi_list),{'V2d','V2v','V3d','V3v','VO1','VO2',...
    'LO1','LO2','TO1','TO2','FOVEAIPS'}),:) = [];

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
    
    [RFs,anti_RFs] = draw_RFs(PRFparams,ROI);
    spec_RFs = RFs(VE>0.1,:);
    
    RF_thresh = 0.001;
    % each stimulus is a gaussian around a polar angle (approx)
    % how aligned each stimulus with each RF of each voxel will be telling
    inRF = (normpdf(0:359,task.stimval,1) * spec_RFs')>RF_thresh;
    outRF = (normpdf(0:359,task.stimval,1) * anti_RFs(VE>0.1,:)')>RF_thresh;

    %initialize saving structure here
    eval(['delay_mean_PSC.' + string(ROI_name) + ' = NaN(max(full.overalltrial),nvoxels);'])

    for sess = 1:sessions
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
            
            % OPTION 1: Divide by the mean of each time course, so the mean over 
            % the 399 TRs
            % PSC = 100 * ((fmriSignal./(nanmean(fmriSignal,2))) - 1);
            % OPTION 2: Divide by the mean of the first few TR's, so all
            % time courses start at 0
            
            data = full(full.run==r&full.sess==sess,:);
            
            % Cycle over trials 
            for t = 1:max(unique(data.trial))
                trial_signal = NaN(nvoxels,55);
                
                %trial_signal(:,1:sum(data.trial==t)) = PSC(:,data.trial==t);
                % Grab TR's inside each trial, in order
                start = fmriSignal(:,data.trial==t); start = start(:,1:4);
                PSC_trial = 100 * ((fmriSignal(:,data.trial==t)./nanmean(start,2)) - 1);
                trial_signal(:,1:sum(data.trial==t)) = PSC_trial;

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
    
    % Get task data
    hard_trials = unique(full.overalltrial(full.cond==2,:));
    hard_idx = sum(trial_index == hard_trials',2)==1;

    figure('Position',[0 0 1200 500])
    subplot(1,2,1)
    % Plot time course of PSC for this ROI, differentiating between
    % conditions
    errorbar(nanmean(all_trial_signal(hard_idx,:)),nanstd(all_trial_signal(hard_idx,:))./sqrt(nvoxels*sess),'Color',condcolors(1,:),'LineWidth',1.25,'DisplayName','Hard trials')
    xticks(0:4:44)
    xticklabels(0:TR*2:TR*22)
    hold on
    errorbar(nanmean(all_trial_signal(~hard_idx,:)),nanstd(all_trial_signal(~hard_idx,:))./sqrt(nvoxels*sess),'Color',condcolors(2,:),'LineWidth',1.25,'DisplayName','Easy trials')
    plot([2.5 2.5],ylim,'k--','LineWidth',1.25,'DisplayName','Stim onset')
    plot([32 32],ylim,'k--','LineWidth',1.25,'DisplayName','Test onset')
    plot([34 34],ylim,'k--','LineWidth',1.25,'DisplayName','Feedback onset')
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 12;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean percent signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    
    subplot(1,2,2)
    eval(['means = delay_mean_PSC.' + string(ROI_name) + ';'])
    % take mean of means, over 2 x 2 in RF vs trial type
    inRF_easy = means(repmat(task.cond==1,1,size(inRF,2))&inRF);
    inRF_hard = means(repmat(task.cond==2,1,size(inRF,2))&inRF);
    outRF_easy = means(repmat(task.cond==1,1,size(outRF,2))&outRF);
    outRF_hard = means(repmat(task.cond==2,1,size(outRF,2))&outRF);
    
    trialtypemeans = [nanmean(means(task.cond==1,:),2) nanmean(means(task.cond==2,:),2)];
    
    errorbar(nanmean(trialtypemeans),[nanstd(trialtypemeans(:,1)) nanstd(trialtypemeans(:,2))]./sqrt(size(means,1)./2), ...
        'ok','LineWidth',2,'DisplayName','Mean over trial types')
    hold on
    errorbar([nanmean(inRF_easy) nanmean(inRF_hard)],[nanstd(inRF_easy) nanstd(inRF_hard)]./sqrt(size(means,1)./2), ...
        'ob','LineWidth',2,'DisplayName','Mean PSC inside RF')
    errorbar([nanmean(outRF_easy) nanmean(outRF_hard)],[nanstd(outRF_easy) nanstd(outRF_hard)]./sqrt(size(means,1)./2), ...
        'or','LineWidth',2,'DisplayName','Mean PSC outside RF')
    
    title(['Mean late delay period activity: ' ROI_name{1}])
    xticks([1 2])
    xticklabels({'Easy trials','Hard trials'})
    ylabel('Mean % signal change')
    ax = gca; ax.FontSize = 12; xlim([0.5 2.5]);
    legend('location','best')
    
    pause(5); saveas(fig,['PSC/' + string(ROI_name{1}) + '.jpg']);
    close all
    
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


