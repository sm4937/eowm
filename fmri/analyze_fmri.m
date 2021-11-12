%% Analyze fMRI data from effort of working memory project
% (EOWM)
% Written by Sarah Master with help from Yuna Kwak, Grace Hallenbeck
% Start date: 2021

clear all

addpath(genpath('/Users/sarah/Documents/MATLAB/fmriTools'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/preproc_mFiles'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/vistasoft_ts'))
% pull up some pre-built functions

condcolors = [190 0 110; 0 110 190]./255;
TR = 0.75;
% Load data up from nifti files

%specify subject ID as string
load('subjinits.mat')
subjlist = [4 5 6];
nsessions = [2 2 2];

subjnum = 6;
sessions = 2;
%specify n sessions for your subject
subject = subjinits{subjnum};

% make sure data accessible
savepath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];
newdatadir = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum) '/fmri/'];
datapath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/'];
mkdir(newdatadir);
subjfiles = ls(newdatadir);
ROIpath = [datapath 'ROIs'];

% if data not already pulled, pull it
for sess = 1:sessions
    if ~contains(subjfiles,['sess' num2str(sess)]) %fmri data has not already been loaded
        pull_fmri_data(subjnum,sess)
    end
end

% Match task data to TR's for all subjects
full = [];
for sii = 1:length(subjlist)
    subjnum = subjlist(sii);
    sessions = nsessions(sii);
    temp = get_full_table(subjnum,sessions);
    full = [full; temp];
end

% Note: To find PRF params for each subject, look for file names
% resembling:
% PRFparams = niftiread('/System/Volumes/Data/d/DATA/data/vRF_tcs/CC/RF1/CC_RF1_vista/RF_ss5_25mm-fFit.nii.gz');


%% Analyze functional data in terms of percent signal change,
% across all PRF-defined ROIs 

% pull an example subject directory, for grabbing a list of ROIs to look at
% across subjects (CC is a good one)
ROIpath = '/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/';
ROI_folder = dir([ROIpath 'CC/ROIs/']);
ROI_list = char(ROI_folder.name); %list all possible ROIs
ROI_list(1:2,:) = []; % get rid of blank indices in folder
%ROI_list = string(ROI_list);

hemi = 'bilat';
% pick the hemisphere you want to examine
hemi_list = ROI_list(contains(string(ROI_list),hemi),:);
% trim ROIs which not all subjects have, for comparison's sake
hemi_list(contains(string(hemi_list),{'V2d','V2v','V3d','V3v','VO1','VO2',...
    'LO1','LO2','TO1','TO2','FOVEAIPS','hV4'}),:) = [];
temp_list(1:4,:) = hemi_list(5:8,:);
temp_list(5:8,:) = hemi_list(1:4,:);
temp_list(9:10,:) = hemi_list(9:10,:); % sort to align by visual hierarchy
hemi_list = temp_list;

%Pre-define some variables that you want averages over
%delay_mean_PSC = struct();

for area = 1:size(hemi_list,1)
    ROI_filename = strrep(hemi_list(area,:),' ','');
    pat = '.(\w*).nii.gz';
    ROI_name = regexp(ROI_filename, pat, 'tokens');
    
    trialtypemeans = []; trialtypeSEMs = [];

    for sii = 1:length(subjlist)
        subjnum = subjlist(sii);
        
        %Initialize these variables
        trial_index = []; all_trial_signal = [];

        % Get task data
        load(['task_subj' num2str(subjnum) '.mat'])
        eval(['task = task_subj' num2str(subjnum) ';'])
        load(['PRFparams_subj' num2str(subjnum) '.mat'])
        
        ROI = niftiRead([ROIpath subjinits{subjnum} '/ROIs/' ROI_filename]);

        % get variance explained from PRF fitting
        VE = PRFparams(:,:,:,2); VE = VE(ROI.data>0);
        selection_cutoff = 0.1;

        [RFs,anti_RFs] = draw_RFs(PRFparams,ROI,'box');
        RF_thresh = 0.001;
        
        inRF = (normpdf(0:359,task.stimval,1) * RFs(VE>=selection_cutoff,:)')>RF_thresh;
        outRF = (normpdf(0:359,task.stimval,1) * anti_RFs(VE>=selection_cutoff,:)')>RF_thresh;
        
        nvoxels_initial = sum(VE>=selection_cutoff);
        %initialize saving structure here
        eval(['delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}) ' = NaN(height(task),nvoxels_initial);'])
        
        sessions = nsessions(sii);
        for sess = 1:sessions
            load(['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum) '/fmri/sess' num2str(sess) '.mat']);
            % load up functional data for that session
            runs = length(funcdata);

            for r = 1:runs %cycle over runs within that session

                all_voxel_tcs = niftiExtract(funcdata{r},ROI);
                % Tommy wrote this function to extract ROI-based functional
                % maps

                all_voxel_tcs = all_voxel_tcs';
                % reshape the voxels into a 2D matrix,
                % where time is the 2nd dimension and voxel # is the first 
                spec_voxel_tcs = all_voxel_tcs(VE>=selection_cutoff,:);
%                 empty = sum(spec_voxel_tcs>0,2)==0;
                
%                 spec_voxel_tcs(empty,:) = []; % this is a patch on a different ROI problem, just to get other things running for now
%                 inRF(:,empty) = []; outRF(:,empty) = [];
                
                nvoxels = size(spec_voxel_tcs,1);

                fmriSignal = spec_voxel_tcs;
                % OPTION 1: Divide by the mean of each time course, so the mean over 
                % the 399 TRs
                % PSC = 100 * ((fmriSignal./(nanmean(fmriSignal,2))) - 1);
                % OPTION 2: Divide by the mean of the first few TR's, so all
                % time courses start at 0 (see below, approx line 150)

                data = full(full.run==r&full.sess==sess&full.subj==subjnum,:);

                % Cycle over trials 
                for t = 1:max(unique(data.trial))
                    trial_signal = NaN(nvoxels,55);

                    % trial_signal(:,1:sum(data.trial==t)) = PSC(:,data.trial==t);
                    
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
                    
                    temp = nanmean(trial_signal(:,late_delay),2)';
                    eval(['delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}{1}) '(overalltrial,1:nvoxels) = temp;'])
                    % save mean PSC in delay period on this trial, for that ROI
                    % resulting struct has trial (rows) x voxel (columns)
                    % measurements of mean delay period PSC
                    
                end % of cycling over trial

            end % of cycling over runs

        end % of cycling over sessions
        
%         eval(['delay_mean_PSC = delay_mean_PSC_subj' num2str(subjnum) ';'])
%         save(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected.mat'],'delay_mean_PSC')
    
        eval(['means = delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}{1}) ';'])
        % take mean of means, over 2 x 2 in RF vs trial type
        inRF_easy = means(repmat(task.cond==1,1,size(inRF,2))&inRF);
        inRF_hard = means(repmat(task.cond==2,1,size(inRF,2))&inRF);
        outRF_easy = means(repmat(task.cond==1,1,size(outRF,2))&outRF);
        outRF_hard = means(repmat(task.cond==2,1,size(outRF,2))&outRF);

        trialtypemeans = [trialtypemeans; nanmean(nanmean(means(task.cond==1,:))) nanmean(nanmean(means(task.cond==2,:))) ...
            nanmean(inRF_easy) nanmean(inRF_hard) nanmean(outRF_easy) nanmean(outRF_hard)];
        trialtypeSEMs = [trialtypeSEMs; nanstd(nanmean(means(task.cond==1,:),2)) nanstd(nanmean(means(task.cond==2,:),2)) ...
            nanstd(inRF_easy) nanstd(inRF_hard) nanstd(outRF_easy) nanstd(outRF_hard)]./sqrt(120);

        hard_trials = task.overalltrial(task.cond==2,:);
        hard_idx = sum(trial_index == hard_trials',2)==1;
    
        % Grab time course of PSC for this ROI, differentiating between
        % conditions, for one subject at a time
        SEMs_hard(sii,:) = nanstd(all_trial_signal(hard_idx,:),1)./sqrt(sum(hard_idx));
        SEMs_easy(sii,:) = nanstd(all_trial_signal(~hard_idx,:),1)./sqrt(sum(~hard_idx));
        tc_hard(sii,:) = nanmean(all_trial_signal(hard_idx,:),1);
        tc_easy(sii,:) = nanmean(all_trial_signal(~hard_idx,:),1);
        
    end % of cycling over participants

    %figure('Position',[0 0 1200 500])
    figure(1)
    subplot(5,2,area)
    % Plot time course of PSC for this ROI, differentiating between
    % conditions
    avg_SEM_hard = nanmean(SEMs_hard(:,1:20));
    avg_SEM_easy = nanmean(SEMs_easy(:,1:20));
    avg_tc_hard = nanmean(tc_hard(:,1:20)); avg_tc_easy = nanmean(tc_easy(:,1:20));

    errorbar(avg_tc_hard,avg_SEM_hard,'Color',condcolors(1,:),'LineWidth',1.25,'DisplayName','Hard trials')
    xticks(0:2:37)
    xticklabels(0:TR*2:TR*37)
    hold on
    errorbar(avg_tc_easy,avg_SEM_easy,'Color',condcolors(2,:),'LineWidth',1.25,'DisplayName','Easy trials')
    %delay period starts 1.5 seconds in
    plot([3 3],ylim,'k--','LineWidth',1.25,'DisplayName','Delay onset')
    %delay period ends at 13.5
    plot([18 18],ylim,'k--','LineWidth',1.25,'DisplayName','Delay end')
    % by 15 seconds in, there's just he ITI, so that should get trimmed out
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 12;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean percent signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    pause(5); saveas(fig,['PSC/Mean PSC all ROIs.jpg']);
    
    figure(2)
    subplot(5,2,area)
    errorbar(nanmean(trialtypemeans(:,1:2),1),[nanmean(trialtypeSEMs(:,1:2),1)], ...
        'ok','LineWidth',2,'DisplayName','Mean over trial types')
    % plotting trial type means & errorbars using averages, average SEMs
    hold on
    errorbar(nanmean(trialtypemeans(:,3:4),1),[nanmean(trialtypeSEMs(:,3:4),1)], ...
        'ob','LineWidth',2,'DisplayName','Mean PSC inside RF')
    errorbar(nanmean(trialtypemeans(:,5:6),1),[nanmean(trialtypeSEMs(:,5:6),1)], ...
        'or','LineWidth',2,'DisplayName','Mean PSC outside RF')
    
    title(['Mean late delay period activity: ' ROI_name{1}])
    xticks([1 2])
    xticklabels({'Easy trials','Hard trials'})
    ylabel('Mean % signal change')
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 14; xlim([0.5 2.5]);
    legend('location','best')
    
    figlabel = ['Subject ' num2str(subjnum)];
    if length(subjlist) > 1; figlabel = ['N = ' num2str(length(subjlist))]; end
    pause(5); saveas(fig,'PSC/Condition differences PSC all ROIs.jpg');
    %saveas(fig,['PSC/' figlabel ':' char(ROI_name{1}) '.jpg']);
        
end % of cycling over areas of interest (ROIs)

%% Make PSC graphs based on RF center distance from stim
clear all
load('PSC/late_delay_mean_PSC_subj4_VEselected.mat')
hemi_list = fieldnames(delay_mean_PSC);
selection_cutoff = 0.1;
condcolors = [190 0 110; 0 110 190]./255;

load('subjinits.mat')
%subjlist = [4 5 6];
subjlist = 6;

for rii = 1:size(hemi_list,1)
    ROI_name = hemi_list{rii};
    ROI_filename = ['bilat.' ROI_name '.nii.gz'];
    
    for sii = 1:length(subjlist)
        subjnum = subjlist(sii);
        ROI = niftiRead(['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subjinits{subjnum} '/ROIs/' ROI_filename]);
        
        load(['PRFparams_subj' num2str(subjnum) '.mat'])
        VE = PRFparams(:,:,:,2);
        
%         centers_x = PRFparams(:,:,:,6); 
%         centers_y = PRFparams(:,:,:,7);
%         % Centers in x and y coordinates
%         ROI_centers = [centers_x(ROI.data>0 & VE >= selection_cutoff) ...
%             centers_y(ROI.data>0 & VE >= selection_cutoff)];
        
        centers = rad2deg(PRFparams(:,:,:,1));
        %centers in degrees
        ROI_centers = round(centers(ROI.data>0&VE>=selection_cutoff,:));
        
        load(['task_subj' num2str(subjnum) '.mat'])
        eval(['task = task_subj' num2str(subjnum) ';'])
        %stimval = [task.x_target task.y_target]; 
        stimval = task.stimval;
        stimval = stimval(sum(~isnan(stimval),2)>0,:);
        cond = task.cond(task.cond>0,:);
        
        load(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected.mat'])
        eval(['means = delay_mean_PSC.' ROI_name ';'])
        
        PSC_sorted_means = []; PSC_binned_means = [];
        for tii = 1:length(stimval)
            
            distances = get_angular_distance(stimval(tii,:),ROI_centers);
            tosort = [means(tii,:)' distances]; 
            
            nbins = 10;
            [y,edges] = histcounts(distances,nbins);
            for bin = 1:nbins
                start_bin = edges(bin);
                end_bin = edges(bin+1);
                PSC_binned_means{bin}(tii,:) = nanmean(tosort((tosort(:,2)<end_bin & tosort(:,2)>=start_bin),1));
            end
            
            % could just sort by absolute distance, or could bin in
            % relation to each other
            sorted_means = sortrows(tosort,2);
            PSC_sorted_means(tii,:) = sorted_means(:,1)';
            % one row per trial, one column per voxel
            
        end % of looping over trials
                
%          hard_means(sii,:) = nanmean(PSC_sorted_means(cond==2,:),1);
%          easy_means(sii,:) = nanmean(PSC_sorted_means(cond==1,:),1);
%          hard_SEMs(sii,:) = nanstd(PSC_sorted_means(cond==2,:),[],1)./sqrt(sum(cond==2));
%          easy_SEMs(sii,:) = nanstd(PSC_sorted_means(cond==1,:),[],1)./sqrt(sum(cond==1));
         
        for bin = 1:nbins
            hard_line(sii,bin) = nanmean(PSC_binned_means{bin}(cond==2,:));
            hard_SEM(sii,bin) = nanstd(PSC_binned_means{bin}(cond==2,:))./sqrt(length(PSC_binned_means{bin}(cond==2,:)));
            easy_line(sii,bin) = nanmean(PSC_binned_means{bin}(cond==1,:));
            easy_SEM(sii,bin) = nanstd(PSC_binned_means{bin}(cond==1,:))./sqrt(length(PSC_binned_means{bin}(cond==1,:)));
        end

         
    end % of looping over subjects

    figure(3)
    subplot(5,2,rii)
    errorbar(nanmean(hard_line,1),nanmean(hard_SEM,1),'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trials')
    hold on
    errorbar(nanmean(easy_line,1),nanmean(easy_SEM,1),'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trials')
    xlabel('Distance of RF center from target stimulus (AU)')
    ylabel('Mean % signal change: Late delay period')
    title([ROI_name])
    legend('Location','Best')
    fig = gcf; fig.Color = 'w';
    ax = gca; ax.FontSize = 12;

    
end % of looping over ROIs

%% Run TAFKAP on data, get estimation accuracy out

% Example subject
subjnum = 6;
% % Grab task data
load(['task_subj' num2str(subjnum) '.mat'])
eval(['task = task_subj' num2str(subjnum) ';'])

%Specify which file you're interested in
spec = '_late_delay';

if contains(ls,['TAFKAP_subj' num2str(subjnum) spec '.mat'])
    % if TAFKAP has been run on this subject already, load that output
    load(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')
    
else %TAFKAP has not been run on this subject yet
    % load in BOLD data (already meaned)
    load(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected.mat'])

    BOLD = delay_mean_PSC;
    TAFKAP_output = run_TAFKAP_eowm(BOLD,task);
    save(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')
    
end % of deciding to run TAFKAP or not


ROI_list = fieldnames(TAFKAP_output);

figure
% Cycle through each ROI and plot estimation accuracy
stimval = task.stimval(~isnan(task.stimval));
for rii = 1:length(ROI_list)
    
    ROI_name = ROI_list{rii};
    
    eval(['estimates = (TAFKAP_output.' ROI_name '.est).*2;'])
    
    subplot(3,4,rii)
    scatter(stimval(1:length(estimates)),estimates,'Filled')
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 12;
    xlabel('Real stimuli'); ylabel('Decoded stimuli')
    title(['Decoding accuracy of ' ROI_name])
    
end

% plot example distributions for easy and hard conditions
task_trimmed = task(~isnan(task.stimval),:);
error_V1 = get_angular_distance((TAFKAP_output.V1.est).*2,task_trimmed.stimval);
error_IPS0 = get_angular_distance((TAFKAP_output.IPS0.est).*2,task_trimmed.stimval);

unc_V1 = TAFKAP_output.V1.unc;
unc_IPS0 = TAFKAP_output.IPS0.unc;

% example_hard_trial = 75;
% example_easy_trial = 164;
example_hard_trial = randsample(find(task_trimmed.cond==2),1);
example_easy_trial = randsample(find(task_trimmed.cond==1),1);
xs = -179:180;

figure; subplot(1,2,1)
V1_hard_dist = normpdf(xs,(error_V1(example_hard_trial)),unc_V1(example_hard_trial));
V1_easy_dist = normpdf(xs,(error_V1(example_easy_trial)),unc_V1(example_easy_trial)); 
plot(xs,V1_hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,V1_easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
title(['Example stimulus representations in V1: Subject ' num2str(subjnum)])
xlabel('Representation error')
ax = gca; ax.FontSize = 12;

subplot(1,2,2)
ips0_hard_dist = normpdf(xs,abs(error_IPS0(example_hard_trial)),unc_IPS0(example_hard_trial));
ips0_easy_dist = normpdf(xs,abs(error_IPS0(example_easy_trial)),unc_IPS0(example_easy_trial)); 
plot(xs,ips0_hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,ips0_easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
title(['Example stimulus representations in IPS0: Subject ' num2str(subjnum)])
xlabel('Representation error')
fig = gcf; fig.Color = 'w'; legend('Location','Best')
ax = gca; ax.FontSize = 12; linkaxes



%% Set up GLM for each subject

initialize_GLM_stuff(subjnum)



