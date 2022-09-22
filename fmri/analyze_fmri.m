%% Analyze fMRI data from effort of working memory project
% (EOWM)
% Written by Sarah Master with help from Yuna Kwak, Grace Hallenbeck
% Start date: 2021

clear all

addpath(genpath('/Users/sarah/Documents/MATLAB/fmriTools'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/preproc_mFiles'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/vistasoft_ts'))
addpath('../')
addpath(genpath('Users/sarah/Documents/MATLAB/RMAOV33/'))
% pull up some pre-built functions

condcolors = [190 0 110; 0 110 190]./255;
TR = 0.75;
nTRs = 399;
% Load data up from nifti files

%specify subject ID as string
load('subjinits.mat')
subjlist = [4 5 6 7 8 9 10 11 12 13 14 16]; 
% there is a subj 15, but her data are all bad (eyetracking messed up, PRF
% maps discontinuous, etc. etc.)
nsessions = [2 2 2 2 2 2 2 2 2 2 2 2];
n = length(subjlist);

subjnum = 16;
sessions = 2;
%specify n sessions for your subject
subject = subjinits{subjnum};

% make sure data accessible
savepath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];
newdatadir = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum) '/fmri/'];
datapath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/'];
%datapath = ['/System/Volumes/Data/d/DATC/datc/eowm_SM/old_preprocessing_fmri/' subject '/'];
mkdir(newdatadir);
% copies of many key data files are stored here now:
% ('/System/Volumes/Data/d/DATA/home/sarah/key_eowm_files')
% make sure they don't blow up if my local computer does for some reason
subjfiles = ls(newdatadir);

% if data not already pulled, pull it
for sess = 1:sessions
    if ~contains(subjfiles,['sess' num2str(sess)]) %fmri data has not already been loaded
        pull_fmri_data(subjnum,sess,datapath)
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

save('current_full.mat','full','-v7.3')
% Note: To find PRF params for each subject, look for file names
% resembling:
% PRFparams = niftiread('/System/Volumes/Data/d/DATA/data/vRF_tcs/CC/RF1/CC_RF1_vista/RF_ss5_25mm-fFit.nii.gz');
% if you ever use y param from PRF params, apply -1 to it. it's a flipped
% output

%% Analyze functional data in terms of percent signal change,
% across all PRF-defined ROIs 

% if ~exist(full)
%     load('current_full.mat')
% end

% V1-V2-V3 BUG REMAINS!!! In 4 line plot (figure 4)
% the problematic subject is subject #12
% it's an EASY OUT RF voxel that is the problem, or maybe a few of them

% pull an example subject directory, for grabbing a list of ROIs to look at
% across subjects (CC is a good one)
subject = 'CC';
%ROIpath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/ROIs'];
ROIpath = ['/System/Volumes/Data/d/DATC/datc/eowm_SM/old_preprocessing_fmri/'];
ROI_folder = dir([ROIpath 'CC/ROIs/']);
ROI_list = char(ROI_folder.name); %list all possible ROIs
ROI_list(1:2,:) = []; % get rid of blank indices in folder
%ROI_list = string(ROI_list);

hemi = 'bilat';
% pick the hemisphere you want to examine
hemi_list = ROI_list(contains(string(ROI_list),hemi),:);
% trim ROIs which not all subjects have, for comparison's sake
hemi_list(contains(string(hemi_list),{'V2v','V3v','VO1','VO2',...
    'LO1','LO2','TO1','TO2','FOVEAIPS','hV4'}),:) = [];
hemi_list(contains(string(hemi_list),[hemi '.V2.nii.gz']),:) = [];
hemi_list(contains(string(hemi_list),[hemi '.V3.nii.gz']),:) = [];
% limit analyses to V2d, not V2 or V2v.
temp_list(1:4,:) = hemi_list(5:8,:);
temp_list(5:8,:) = hemi_list(1:4,:);
temp_list(9:10,:) = hemi_list(9:10,:); % sort to align by visual hierarchy
hemi_list = temp_list;
spec = '';

% mega_ROIs = {'V1_V2_V3_V3AB', 'IPS0_IPS1_IPS2_IPS3', 'iPCS_sPCS'}';
% hemi_list = mega_ROIs;
% spec = '_megaROIs';

fovea_ROIs = {'V1_V2d_V3d', 'V3AB', 'IPS0_IPS1', 'IPS2_IPS3', 'iPCS', 'sPCS'}';
hemi_list = fovea_ROIs;
spec = '_foveaROIs';

%Pre-define some variables that you want averages over
%delay_mean_PSC = struct();
SEMs_hard = NaN(length(subjlist),55);
SEMs_easy = NaN(length(subjlist),55);
tc_hard = NaN(length(subjlist),55);
tc_easy = NaN(length(subjlist),55);

for area = 1:size(hemi_list,1)
    clear ROI ROI_name
    if contains(spec,'fovea') % foveal ROI
        ROI_name{1} = fovea_ROIs{area};
        all_names = strsplit(ROI_name{1},'_');
        ROI_filename = [repmat('bilat.',length(all_names),1) + string(all_names') + repmat('.nii.gz',length(all_names),1)];
        %ROI_filename = ROI_filename{1};
    elseif contains(spec,'mega') % mega ROI (built up of many individual ROIs)
        ROI_name{1} = mega_ROIs{area};
        all_names = strsplit(ROI_name{1},'_');
        ROI_filename = [repmat('bilat.',length(all_names),1) + string(all_names') + repmat('.nii.gz',length(all_names),1)];
        ROI_filename = ROI_filename{1};
    else % just one ROI at a time, no combining across PRF-defined areas
        ROI_filename = strrep(hemi_list(area,:),' ','');
        pat = '.(\w*).nii.gz';
        ROI_name = regexp(ROI_filename, pat, 'tokens');
    end
    
    trialtypemeans = NaN(length(subjlist),6); trialtypeSEMs = NaN(length(subjlist),6);

    for sii = 1:length(subjlist)
        subjnum = subjlist(sii);
        
        %Initialize these variables
        trial_index = []; all_trial_signal = []; tc_hard_inRF_subj = [];
        tc_easy_inRF_subj = []; tc_easy_outRF_subj = []; tc_hard_outRF_subj = [];

        % Get task data
        load(['task_subj' num2str(subjnum) '.mat'])
        eval(['task = task_subj' num2str(subjnum) ';'])
        load(['PRFparams_subj' num2str(subjnum) '.mat'])
        
        ROI_path_subj = [ROIpath subjinits{subjnum} '/ROIs/'];
        
        %initialize saving structure here
        nvoxels_initial = 1;
        % this variable should be dictated by subj CC
        eval(['delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}) ' = NaN(height(task),nvoxels_initial);'])
        
        has_all_ROIs = true;
        for jj = 1:length(ROI_filename)
            has_all_ROIs = has_all_ROIs & contains(ls(ROI_path_subj),ROI_filename(jj,:));
        end
        if subjnum==12 && area==1 %patch this bug this way
            % clay thinks I've grabbed a blood vessel by accident in my ROI
            % in this person (it's 2-3 voxels I think)
            has_all_ROIs = false;
        end
        
        if has_all_ROIs
            % exclude subjects from means who don't have bilateral ROI in
            % this area
            ROI = niftiRead([ROI_path_subj ROI_filename{1}]);
            % initialize with one template ROI

            ROI.data = false(size(PRFparams(:,:,:,1)));
            for rii = 1:length(ROI_filename)
                temp = niftiRead([ROIpath subjinits{subjnum} '/ROIs/' ROI_filename{rii}]);
                ROI.data = ROI.data | temp.data;
            end
            % get variance explained from PRF fitting
            VE = PRFparams(:,:,:,2); VE = VE(ROI.data>0);
            selection_cutoff = 0.1;

            [RFs,anti_RFs] = draw_RFs(PRFparams,ROI,'box');
            RF_thresh = 0.001; % consider increasing this threshold

            inRF = (normpdf(0:359,task.stimval,1) * RFs(VE>=selection_cutoff,:)')>RF_thresh;
            outRF = (normpdf(0:359,task.stimval,1) * anti_RFs(VE>=selection_cutoff,:)')>RF_thresh;

            inRF_hard = NaN(1,length(task.stimval));
            inRF_easy = inRF_hard;
            outRF_hard = inRF_hard;
            outRF_easy = inRF_hard;
            
            hard_storage_long = NaN(length(task.stimval),1);
            
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
                    empty = sum(spec_voxel_tcs<100,2)==nTRs;

                    fmriSignal = spec_voxel_tcs(~empty,:);
                    
                    inRF = (normpdf(0:359,task.stimval,1) * RFs(VE>=selection_cutoff,:)')>RF_thresh;
                    outRF = (normpdf(0:359,task.stimval,1) * anti_RFs(VE>=selection_cutoff,:)')>RF_thresh;
                    inRF(:,empty) = []; outRF(:,empty) = [];
                    
                    nvoxels = size(fmriSignal,1);
                    
                    % OPTION 1: Divide by the mean of each time course, so the mean over 
                    % the 399 TRs
                    % PSC = 100 * ((fmriSignal./(nanmean(fmriSignal,2))) - 1);
                    % OPTION 2: Divide by the mean of the first few TR's, so all
                    % time courses start at 0 (see below, approx line 150)

                    data = full(full.run==r&full.sess==sess&full.subj==subjnum,:);

                    % Cycle over trials 
                    for t = 1:max(unique(data.trial))
                        trial_signal = NaN(nvoxels,55);

                        overalltrial = unique(data.overalltrial(data.trial==t));
                        trial_index = [trial_index; repmat(overalltrial,size(trial_signal,1),1)];
                        % Keep an index of which trial is in which row (remember that
                        % each row is a voxel's time course inside that trial, and that
                        % we're averaging across voxels, over all trials. In the
                        % future you may want to use this index to compare percent
                        % signal change across trial types.)

                        % Grab TR's inside each trial, in order
                        start = fmriSignal(:,data.trial==t); start = start(:,1:4);
                        PSC_trial = 100 * ((fmriSignal(:,data.trial==t)./nanmean(start,2)) - 1);
                        trial_signal(:,1:sum(data.trial==t)) = PSC_trial;
                        
                        %tracking some variables - why the huge outliers in
                        %subject 12 V1 timecourse?
                        nanmean(PSC_trial(:));
                        beep = sort(start(:),'descend'); beep(end-150:end-50);
                        
                        all_trial_signal = [all_trial_signal; trial_signal];
                        
                        hard = unique(data.cond(data.trial==t,:))==2;
                        hard_storage_long(overalltrial) = hard;
                        if hard
                            tc_hard_inRF_subj = [tc_hard_inRF_subj; trial_signal(inRF(overalltrial,:)',:)];
                            tc_hard_outRF_subj = [tc_hard_outRF_subj; trial_signal(outRF(overalltrial,:)',:)];
                            inRF_hard(:,overalltrial) = nanmean(nanmean(trial_signal(inRF(overalltrial,:)',:)));
                            outRF_hard(:,overalltrial) = nanmean(nanmean(trial_signal(outRF(overalltrial,:)',:)));
                        else
                            tc_easy_inRF_subj = [tc_easy_inRF_subj; trial_signal(inRF(overalltrial,:)',:)];
                            tc_easy_outRF_subj = [tc_easy_outRF_subj; trial_signal(outRF(overalltrial,:)',:)];
                            inRF_easy(:,overalltrial) = nanmean(nanmean(trial_signal(inRF(overalltrial,:)',:)));
                            outRF_easy(:,overalltrial) = nanmean(nanmean(trial_signal(outRF(overalltrial,:)',:)));
                        end
                        
                        % Grab percent signal change in that voxel over each of those
                        % trials

                        delay = (data.epoch(data.trial==t)==3)';
                        late_delay = data.delay_times(data.trial==t) >= 6;

                        temp = nanmean(trial_signal(:,late_delay),2)';
                        eval(['delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}) '(overalltrial,1:nvoxels) = temp;'])
                        % save mean PSC in delay period on this trial, for that ROI
                        % resulting struct has trial (rows) x voxel (columns)
                        % measurements of mean delay period PSC

                    end % of cycling over trial

                end % of cycling over runs

            end % of cycling over sessions

            eval(['delay_mean_PSC = delay_mean_PSC_subj' num2str(subjnum) ';'])
            save(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected' spec '.mat'],'delay_mean_PSC')
            eval(['means = delay_mean_PSC_subj' num2str(subjnum) '.' char(ROI_name{1}) ';'])
            % take mean of means, over 2 x 2 in RF vs trial type

             trialtypemeans(sii,:) = [nanmean(nanmean(means(task.cond==1,:))) nanmean(nanmean(means(task.cond==2,:))) ...
                 nanmean(inRF_easy) nanmean(inRF_hard) nanmean(outRF_easy) nanmean(outRF_hard)];
             trialtypeSEMs(sii,:) = [nanstd(nanmean(means(task.cond==1,:),2)) nanstd(nanmean(means(task.cond==2,:),2)) ...
                 nanstd(inRF_easy) nanstd(inRF_hard) nanstd(outRF_easy) nanstd(outRF_hard)]...
                 ./sqrt([sum(hard_storage_long==0) sum(hard_storage_long==1) sum(hard_storage_long==0) sum(hard_storage_long==1) sum(hard_storage_long==0) sum(hard_storage_long==1)]);


            hard_trials = task.overalltrial(task.cond==2,:);
            hard_idx = sum(trial_index == hard_trials',2)==1;

            % Grab time course of PSC for this ROI, differentiating between
            % conditions, for one subject at a time
            SEMs_hard(sii,:) = nanstd(all_trial_signal(hard_idx,:),1)./sqrt(sum(hard_storage_long==1));
            SEMs_easy(sii,:) = nanstd(all_trial_signal(~hard_idx,:),1)./sqrt(sum(hard_storage_long==0));
            tc_hard(sii,:) = nanmean(all_trial_signal(hard_idx,:),1);
            tc_easy(sii,:) = nanmean(all_trial_signal(~hard_idx,:),1);
            
            % Now grab time course of hard/easy trials in/out receptive
            % field voxels
            tc_hard_inRF(sii,:) = nanmean(tc_hard_inRF_subj);
            tc_easy_inRF(sii,:) = nanmean(tc_easy_inRF_subj);
            tc_hard_outRF(sii,:) = nanmean(tc_hard_outRF_subj);
            tc_easy_outRF(sii,:) = nanmean(tc_easy_outRF_subj);
            
            eval(['delay_trialtype_PSC.' char(ROI_name{1}) ' = trialtypemeans(sii,:);'])
            eval(['delay_trialtype_PSC.' char(ROI_name{1}) '_SEM = trialtypeSEMs(sii,:);'])
            save(['PSC/late_delay_mean_PSC_trialtypes_subj' num2str(subjnum) '_VEselected_foveaROIs.mat'],'delay_trialtype_PSC')

        end % end of cutting subjects out who don't have this ROI
        
    end % of cycling over participants (sii loop)

    %figure('Position',[0 0 1200 500])
    % Plot time course of PSC for this ROI, differentiating between
    % HARD and EASY trials
    figure(1)
    subplot(3,2,area)
    % Plot time course of PSC for this ROI, differentiating between
    % conditions
    % collapsed across in/out RF
%     avg_tc_hard = nanmean(tc_hard(:,1:20)); avg_tc_easy = nanmean(tc_easy(:,1:20));
%     SEM_hard = nanstd(tc_hard(:,1:20))/sqrt(n);
%     SEM_easy = nanstd(tc_easy(:,1:20))/sqrt(n);
    % in RF only
    avg_tc_hard = nanmean(tc_hard_inRF(:,1:20)); avg_tc_easy = nanmean(tc_easy_inRF(:,1:20));
    SEM_hard = nanstd(tc_hard_inRF(:,1:20))/sqrt(n);
    SEM_easy = nanstd(tc_easy_inRF(:,1:20))/sqrt(n);
    
    hold on
    btwn_fill = [avg_tc_hard + SEM_hard, fliplr(avg_tc_hard - SEM_hard)];
    fill_xs = [1:20, fliplr(1:20)];
    fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials');
    btwn_fill = [avg_tc_easy + SEM_easy, fliplr(avg_tc_easy-SEM_easy)];
    fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials');
    plot(1:20,avg_tc_easy,'Color',condcolors(2,:),'LineWidth',2)
    plot(1:20,avg_tc_hard,'Color',condcolors(1,:),'LineWidth',2)

    %errorbar(avg_tc_hard,SEM_hard,'Color',condcolors(1,:),'LineWidth',1.25,'DisplayName','Hard trials')
    %errorbar(avg_tc_easy,SEM_easy,'Color',condcolors(2,:),'LineWidth',1.25,'DisplayName','Easy trials')
    xticks(0:2:37)
    xticklabels(0:TR*2:TR*37)
    %delay period starts 1.5 seconds in
    plot([3 3],ylim,'k--','LineWidth',1.25,'DisplayName','Delay onset')
    %delay period ends at 13.5
    plot([18 18],ylim,'k--','LineWidth',1.25,'DisplayName','Delay end')
    % by 15 seconds in, there's just he ITI, so that should get trimmed out
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 14;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean % signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    if size(hemi_list,1) == area
        % last ROI has been plotted
        saveas(fig,['PSC/Mean PSC in RF only.jpg']);
    end
    pause(1)
    
    figure(2)
    subplot(3,2,area)
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
    if size(hemi_list,1) == area
        saveas(fig,'PSC/Condition differences PSC.jpg');
    end
    pause(1)
    

    figure(3)
    subplot(3,2,area)
    % Plot time course of PSC for this ROI, differentiating between
    % IN and OUT RF voxels
    neutral_color = [50 50 50]./255;
    avg_tc_in = nanmean(tc_hard_inRF(:,1:20)); avg_tc_out = nanmean(tc_hard_outRF(:,1:20));
    SEM_in = nanstd(tc_hard_inRF(:,1:20))/sqrt(n);
    SEM_out = nanstd(tc_hard_outRF(:,1:20))/sqrt(n);
    
    hold on
    btwn_fill = [avg_tc_in + SEM_in, fliplr(avg_tc_in - SEM_in)];
    fill_xs = [1:20, fliplr(1:20)];
    fill(fill_xs,btwn_fill,neutral_color,'linestyle','none','facealpha',0.2,'DisplayName','In RF');
    btwn_fill = [avg_tc_out + SEM_out, fliplr(avg_tc_out-SEM_out)];
    fill(fill_xs,btwn_fill,neutral_color,'linestyle','none','facealpha',0.2,'DisplayName','Out of RF');
    plot(1:20,avg_tc_in,'Color',neutral_color,'LineWidth',2)
    plot(1:20,avg_tc_out,'--','Color',neutral_color,'LineWidth',2)
    xticks(0:2:37)
    xticklabels(0:TR*2:TR*37)
    %delay period starts 1.5 seconds in
    plot([3 3],ylim,'k--','LineWidth',1.25,'DisplayName','Delay onset')
    %delay period ends at 13.5
    plot([18 18],ylim,'k--','LineWidth',1.25,'DisplayName','Delay end')
    % by 15 seconds in, there's just he ITI, so that should get trimmed out
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 14;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean % signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    if size(hemi_list,1) == area
        % last ROI has been plotted
        saveas(fig,['PSC/Mean PSC in vs out RF (hard only).jpg']);
    end
    pause(1)
    
    % superimpose ALL conditions here
    figure(4)
    subplot(3,2,area)
    % Plot time course of PSC for this ROI, differentiating between
    % conditions
    SEM_hard_inRF = nanstd(tc_hard_inRF(:,1:20))/sqrt(n);
    SEM_easy_inRF = nanstd(tc_easy_inRF(:,1:20))/sqrt(n);
    SEM_hard_outRF = nanstd(tc_hard_outRF(:,1:20))/sqrt(n);
    SEM_easy_outRF = nanstd(tc_easy_outRF(:,1:20))/sqrt(n);
    
    avg_tc_hard_inRF = nanmean(tc_hard_inRF(:,1:20)); 
    avg_tc_easy_inRF = nanmean(tc_easy_inRF(:,1:20));
    avg_tc_hard_outRF = nanmean(tc_hard_outRF(:,1:20)); 
    avg_tc_easy_outRF = nanmean(tc_easy_outRF(:,1:20));
    
    hold on
    btwn_fill = [avg_tc_hard_inRF + SEM_hard_inRF, fliplr(avg_tc_hard_inRF-SEM_hard_inRF)];
    fill_xs = [1:20, fliplr(1:20)];
    fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials, in RF');
    xticks(0:2:37)
    xticklabels(0:TR*2:TR*37)
    
    btwn_fill = [avg_tc_hard_outRF + SEM_hard_outRF, fliplr(avg_tc_hard_outRF-SEM_hard_outRF)];
    fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2,'DisplayName','Hard trials, out RF');
    btwn_fill = [avg_tc_easy_inRF + SEM_easy_inRF, fliplr(avg_tc_easy_inRF-SEM_easy_inRF)];
    fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials, in RF');
    btwn_fill = [avg_tc_easy_outRF + SEM_easy_outRF, fliplr(avg_tc_easy_outRF-SEM_easy_outRF)];
    fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2,'DisplayName','Easy trials, out RF');
    plot(1:20,avg_tc_hard_outRF,'--','Color',condcolors(1,:),'LineWidth',2)
    plot(1:20,avg_tc_easy_outRF,'--','Color',condcolors(2,:),'LineWidth',2)
    plot(1:20,avg_tc_easy_inRF,'Color',condcolors(2,:),'LineWidth',2)
    plot(1:20,avg_tc_hard_inRF,'Color',condcolors(1,:),'LineWidth',2)

    %errorbar(avg_tc_hard_inRF,SEM_hard_inRF,'Color',condcolors(1,:),'LineWidth',3,'DisplayName','Hard trials, in RF')
    %errorbar(avg_tc_easy_inRF,SEM_easy_inRF,'Color',condcolors(2,:),'LineWidth',3,'DisplayName','Easy trials, in RF')
    %errorbar(avg_tc_hard_outRF,SEM_hard_outRF,'--','Color',condcolors(1,:),'LineWidth',3,'DisplayName','Hard trials, out RF')
    %errorbar(avg_tc_easy_outRF,SEM_easy_outRF,'--','Color',condcolors(2,:),'LineWidth',3,'DisplayName','Easy trials, out RF')
    %delay period starts 1.5 seconds in
    plot([3 3],ylim,'k--','LineWidth',3,'DisplayName','Delay onset')
    %delay period ends at 13.5
    plot([18 18],ylim,'k--','LineWidth',3,'DisplayName','Delay end')
    % by 15 seconds in, there's just he ITI, so that should get trimmed out
    fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 14;
    title(['Mean Signal Change over Trials: ' ROI_name{1}])
    ylabel('Mean % signal change')
    xlabel('Time (seconds)')
    legend('boxoff')
    if size(hemi_list,1) == area
        saveas(fig,['PSC/Mean PSC in and out RF timecourse.jpg']);
    end
    pause(1)
    
        
end % of cycling over areas of interest (ROIs)


%% Run statistics on fMRI timecourses - condition differences

% take mean delay period activity in:
% one number one subject one ROI (12 x 6 per cell)
% - in RF hard 
% - in RF easy
% - out RF hard
% - out RF easy
% - times 6 for each ROI

    errorbar(nanmean(trialtypemeans(:,1:2),1),[nanmean(trialtypeSEMs(:,1:2),1)], ...
        'ok','LineWidth',2,'DisplayName','Mean over trial types')
    % plotting trial type means & errorbars using averages, average SEMs
    hold on
    errorbar(nanmean(trialtypemeans(:,3:4),1),[nanmean(trialtypeSEMs(:,3:4),1)], ...
        'ob','LineWidth',2,'DisplayName','Mean PSC inside RF')
    errorbar(nanmean(trialtypemeans(:,5:6),1),[nanmean(trialtypeSEMs(:,5:6),1)], ...
        'or','LineWidth',2,'DisplayName','Mean PSC outside RF')

% run 3-way ANOVA [within ROI]
% rma0v3
% obtain true F value across all subjects

% now, shuffle WITHIN ROI
% hard & easy trial labels (shuffle of one column in ANOVA table)
% obtain new F value
% rinse, repeat N_permutations times

% obtain final percentage (p-value) of null F values >= true F value
% sum F_boot > F_true
% divide by N_permutations

%% Run TAFKAP on data, get estimation accuracy out

% Example subject
subjnum = 16;
% % Grab task data
load(['task_subj' num2str(subjnum) '.mat'])
eval(['task = task_subj' num2str(subjnum) ';'])
condcolors = [190 0 110; 0 110 190]./255;
neutral_color = [50 50 50]./255;


%Specify which file you're interested in
%spec = '_megaROIs';
%spec = '_late_delay';
spec = '_foveaROIs';

if contains(ls,['TAFKAP_subj' num2str(subjnum) spec '.mat'])
    % if TAFKAP has been run on this subject already, load that output
    load(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')
    
else %TAFKAP has not been run on this subject yet
    % load in BOLD data (already meaned)
    load(['PSC/late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected' spec '.mat'])

    BOLD = delay_mean_PSC;
    TAFKAP_output = run_TAFKAP_eowm(BOLD,task);
    save(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')
    
end % of deciding to run TAFKAP or not

% % % MAKE TAFKAP PLOTS % % %
ROI_list = fieldnames(TAFKAP_output);
%ROI_list = {'V1_V2d_V3d','V3AB','IPS0_IPS1','IPS2_IPS3','iPCS','sPCS'};
ROI_labels = ROI_list;
for ii = 1:length(ROI_labels)
    label = strrep(ROI_labels{ii},'_','-');
    ROI_labels{ii} = label;
end
%ROI_labels{1} = 'V1-V2-V3';

% PLOT DECODING ACCURACY FOR EVERY ROI %
% Cycle through each ROI and plot estimation accuracy
stimval = task.stimval(~isnan(task.stimval));

figure; 
plotnum = [1 2 3 7 8 9];
for rii = 1:length(ROI_list)
    
    ROI_name = ROI_list{rii};
    
    eval(['estimates = (TAFKAP_output.' ROI_name '.est).*2;'])
    
    subplot(4,3,plotnum(rii))
    scatter(stimval(1:length(estimates)),estimates,'Filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.2)
    ylabel('Decoded stimuli'); xlabel('Real stimuli')
    title(ROI_labels{rii})
    [fig,ax] = clean_fig();
    
    subplot(4,3,plotnum(rii)+3)
    distances = get_angular_distance(estimates,stimval(1:length(estimates)));
    histogram(distances)
    xlabel('Absolute decoding error')
    [fig,ax] = clean_fig();
    
    if rii > 6
        % last ROI in the list
        xlabel('Real stimuli')
    end
    
end

%% Plot example distributions for easy and hard conditions
% Averaged across available subjects


subjlist = [4 5 6 7 8 9 10 11 12 13 14 15];
n = length(subjlist);
task_trimmed = task(~isnan(task.stimval),:);
% example_hard_trial = 75;
% example_easy_trial = 164;
example_hard_trial = randsample(find(task_trimmed.cond==2),1);
example_easy_trial = randsample(find(task_trimmed.cond==1),1);
xs = -179:180;

error = NaN(240,n);
unc = NaN(240,n);

% PLOT GROUP MEAN DECODING ACCURACY FOR EVERY ROI & EXAMPLE POSTERIOR DISTRIBUTIONS %

figure; count = 0; xtick_list = [];
for ROI = 1:length(ROI_list)
    
    ROI_name = ROI_list{ROI};
    
    for sii = 1:n
        
        subjnum = subjlist(sii);
        load(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')

        load(['task_subj' num2str(subjnum) '.mat'])
        eval(['task = task_subj' num2str(subjnum) ';'])
        task_trimmed = task(~isnan(task.stimval),:);
        
        eval(['error(1:length(task_trimmed.stimval),sii) = get_angular_distance((TAFKAP_output.' ROI_name '.est).*2,task_trimmed.stimval);'])
        eval(['unc(1:length(task_trimmed.stimval),sii) = TAFKAP_output.' ROI_name '.unc;'])
        
        mean_errors(sii,1) = nanmean(error(task_trimmed.cond==1,sii));
        mean_errors(sii,2) = nanmean(error(task_trimmed.cond==2,sii));
        mean_uncs(sii,1) = nanmean(unc(task_trimmed.cond==1,sii));
        mean_uncs(sii,2) = nanmean(unc(task_trimmed.cond==2,sii));
        
        eval(['est = TAFKAP_output.' ROI_name '.est;'])
        ROI_fidelity(ROI,sii) = corr(task_trimmed.stimval(1:length(est)),est);
        param_tradeoff(ROI,sii) = corr(error(:,sii),unc(:,sii));
        

    end

    if strcmp(ROI_name,'V3AB')
        V3AB_error = mean_errors;
        V3AB_unc = mean_uncs;
    end
    
    xtick_list = [xtick_list; count+1.5];
    
    count = count + 1;
    subplot(3,1,1)
    errorbar(count,nanmean(mean_errors(:,1)),nanstd(mean_errors(:,1))./sqrt(n),'LineWidth',1.5,'Color',condcolors(2,:))
    hold on
    plot(count,nanmean(mean_errors(:,1)),'o','LineWidth',2,'Color',condcolors(2,:))
    errorbar(count+1,nanmean(mean_errors(:,2)),nanstd(mean_errors(:,2))./sqrt(n),'LineWidth',1.5,'Color',condcolors(1,:))
    plot(count+1,nanmean(mean_errors(:,2)),'o','LineWidth',2,'Color',condcolors(1,:))
    [h,p] = ttest(mean_errors(:,1),mean_errors(:,2));
    if p < 0.05
        yvalue = max([nanmean(mean_errors(:,1)) nanmean(mean_errors(:,2))])+7;
        scatter(count+0.5,yvalue,'k*')
    end
    title('Decoding error')
    ylabel('Error (degrees)')
    [fig,ax] = clean_fig();
    xticks(xtick_list)
    xticklabels(ROI_labels)
        
    subplot(3,1,2)
    errorbar(count,nanmean(mean_uncs(:,1)),nanstd(mean_uncs(:,1))./sqrt(n),'LineWidth',1.5,'Color',condcolors(2,:))
    hold on
    plot(count,nanmean(mean_uncs(:,1)),'o','LineWidth',2,'Color',condcolors(2,:))
    errorbar(count+1,nanmean(mean_uncs(:,2)),nanstd(mean_uncs(:,2))./sqrt(n),'LineWidth',1.5,'Color',condcolors(1,:))
    plot(count+1,nanmean(mean_uncs(:,2)),'o','LineWidth',2,'Color',condcolors(1,:))
    [h,p] = ttest(mean_uncs(:,1),mean_uncs(:,2));
    if p < 0.05
        yvalue = max([nanmean(mean_uncs(:,1)) nanmean(mean_uncs(:,2))])+7;
        scatter(count+0.5,yvalue,'k*')
    end
    title(['Decoding uncertainty'])
    ylabel('Uncertainty (degrees)')
    [fig,ax] = clean_fig();
    count = count+3;
    xticks(xtick_list)
    xticklabels(ROI_labels)
    

end %of cycling over ROIs

subjnum = 6; ind = subjnum==subjlist;
subplot(3,1,3)
hard_dist = normpdf(xs,V3AB_error(ind,2),V3AB_unc(ind,2));
easy_dist = normpdf(xs,V3AB_error(ind,1),V3AB_unc(ind,1)); 
%hard_dist = normpdf(xs,nanmean(V3AB_error(:,2)),nanmean(V3AB_unc(:,2)));
%easy_dist = normpdf(xs,nanmean(V3AB_error(:,1)), nanmean(V3AB_unc(:,1)));
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
legend('Location','Best')
title(['Example stimulus representations in V3AB'])
xlabel('Representation error (degrees)')
ylabel('Posterior probability')
yticks([])
[fig,ax] = clean_fig();


% Determine fidelity of representations in each ROI 
figure(8); figure(9);
for ROI = 1:length(ROI_list)
    
    figure(8)
    scatter(ROI.*ones(size(ROI_fidelity,2),1),ROI_fidelity(ROI,:),'Filled')
    hold on
    star_height = 0.7;
    [h,p] = ttest(ROI_fidelity(ROI,:));
    if p < 0.05
        plot(ROI,star_height,'k*')
    end
    errorbar(ROI,nanmean(ROI_fidelity(ROI,:)),nanstd(ROI_fidelity(ROI,:))/sqrt(n),'k','LineWidth',1.5)
    plot(xlim,[0 0],'--','Color',neutral_color)
    [fig,ax] = clean_fig();
    title('Decoding performance by ROI')
    ylabel(sprintf('Correlation \n(Decoded location vs real location)'))
    xlim([0.5 6.5])
    ylim([-0.2 0.8])
    xticklabels(ROI_labels)
    
    figure(9)
    scatter(ROI.*ones(size(param_tradeoff,2),1),param_tradeoff(ROI,:),'Filled')
    hold on
    star_height = 0.8;
    [h,p] = ttest(param_tradeoff(ROI,:));
    if p < 0.05
        plot(ROI,star_height,'*')
    end
    errorbar(ROI,nanmean(param_tradeoff(ROI,:)),nanstd(param_tradeoff(ROI,:))/sqrt(n),'k','LineWidth',1.5)
    xlim([0.5 6.5])
    ylim([-0.2 0.8])
   
        
end
[fig,ax] = clean_fig();
title('Are error & unc correlated? (TAFKAP check)')
ylabel('TAFKAP error versus uncertainty (r)')
xticks(1:length(ROI_list))
xticklabels(ROI_labels)
xlim([0.5 length(ROI_list)+0.5])
ylim([-0.2 star_height+0.1])

%% Supplementary analyses of TAFKAP read-outs

% are TAFKAP read-outs biased towards or away from any areas in particular

count = 0; xval_ROI = [-0.25 0.25];
figure(10); fig = gcf; fig.Color = 'w';
figure(11); fig = gcf; fig.Color = 'w';
figure(12); fig = gcf; fig.Color = 'w';
for ROI = 1:length(ROI_list)
    
    ROI_name = ROI_list{ROI};
    meridian_biases_hard = NaN(n,4); meridian_biases_easy = NaN(n,4);
    cw_biases_easy = NaN(n,1); cw_biases_hard = NaN(n,1);
    response_biases = NaN(n,1);

    for sii = 1:n

        subjnum = subjlist(sii);
        load(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')
        
        eval(['has_ROI = sum(~isnan(TAFKAP_output.' ROI_name '.est))>0;'])

        if has_ROI

            load(['task_subj' num2str(subjnum) '.mat'])
            eval(['task = task_subj' num2str(subjnum) ';'])
            task_trimmed = task(~isnan(task.stimval),:);

            eval(['[cw_biases_onesubj,meridian_biases_onesubj] = get_angular_biases(TAFKAP_output.' ROI_name '.est.*2,task_trimmed.stimval);'])

            cw_biases_easy(sii) = nanmean(cw_biases_onesubj(task_trimmed.cond==1,:))-1.5;
            cw_biases_hard(sii) = nanmean(cw_biases_onesubj(task_trimmed.cond==2,:))-1.5;
            cw_biases(sii,:) = nanmean(cw_biases_onesubj)-1.5;
            
            response_biases(sii) = nanmean(task.resp)-1.5;
            %center on 0

            meridian_biases_easy(sii,:) = nanmean(meridian_biases_onesubj(task_trimmed.cond==1,:));
            meridian_biases_hard(sii,:) = nanmean(meridian_biases_onesubj(task_trimmed.cond==2,:));
        
        end % of ROI if statement
        

    end % of subject-level loop

    figure(10)
    hold on
    xval_ROI = xval_ROI+1;
    scatter(xval_ROI(1)*ones(n,1),cw_biases_easy,[],condcolors(2,:),'Filled')
    errorbar(xval_ROI,[nanmean(cw_biases_easy) nanmean(cw_biases_hard)], ...
        [nanstd(cw_biases_easy) nanstd(cw_biases_hard)]/sqrt(n),'k','LineWidth',1.5)
    scatter(xval_ROI(2)*ones(n,1),cw_biases_hard,[],condcolors(1,:),'Filled')
    plot(xlim,[0 0],'k--')
    ax = gca; ax.FontSize = 14;
    %count = count+3;
    xticks([1:length(ROI_labels)])
    xticklabels(ROI_labels)
    ylabel('Clockwise bias')
    hold off
    
    figure(11)
    subplot(3,2,ROI)
    % 6 ROIs to plot scatters for
    scatter(response_biases,cw_biases,'Filled')
    toinclude = ~isnan(response_biases)&~isnan(cw_biases);
    [r,p] = corr(response_biases(toinclude),cw_biases(toinclude));
    % individual to each ROI
    text = ROI_labels{ROI};
    if p < 0.05
        text = [ROI_labels{ROI} ', r = ' num2str(r)];
        lsline
    end
    title(text)
    ylabel('TAFKAP readout bias')
    xlabel('Subject response bias')
    ax = gca; ax.FontSize = 14;
    
    figure(12)
    subplot(3,2,ROI)
    
    xval = [-0.25 0.25];
    hold on
    for ii = 1:4
        xval = xval+1;
        scatter(xval(1)*ones(n,1),meridian_biases_easy(:,ii),[],condcolors(2,:),'Filled')
        errorbar(xval,[nanmean(meridian_biases_easy(:,ii)) nanmean(meridian_biases_hard(:,ii))], ...
            [nanstd(meridian_biases_easy(:,ii)) nanstd(meridian_biases_hard(:,ii))]/sqrt(n),'k','LineWidth',1.5)
        scatter(xval(2)*ones(n,1),meridian_biases_hard(:,ii),[],condcolors(1,:),'Filled')
    end
    plot(xlim,[0 0],'k--')
    title(ROI_name)
    xlabel('Meridian')
    ylabel('Degrees of bias towards meridian')
    xticks([1:4])
    xlim([0.70 4.30])
    xticklabels({'0 degrees','90 degrees','180 degrees','270 degrees'})
    xtickangle(45)
    ax = gca; ax.FontSize = 14;

end



% for each subject, correlate pupil size and memory precision
% memory precision alreayd pulled trial-by-trial above
load('../data/pupil.mat')

tallyerror = 0; tallyunc = 0;

for sii = 1:length(subjlist)
    
    subj = subjlist(sii);
    idx = pupil.subj==subj;
    
    tocorr = [pupil.pres(idx,:)',pupil.earlydelay(idx,:)',pupil.delay(idx,:)',pupil.latedelay(idx,:)',error(:,sii),unc(:,sii)];
    tocorr(sum(isnan(tocorr),2)>0,:) = [];
    % delete rows where pupil size is NaN
    
    [allRs{sii},allPs{sii}] = corr(tocorr);
    
    tallyerror = tallyerror + sum(allPs{sii}(5,1:4)<0.01);
    tallyunc = tallyunc + sum(allPs{sii}(6,1:4)<0.01);
    
end

disp([num2str(tallyerror) ' significant correlations w/ representational error out of ' num2str(sii) ' subjects, 4 pupil measures.']) 
disp([num2str(tallyunc) ' significant correlations w/ representational uncertainty out of ' num2str(sii) ' subjects, 4 pupil measures.']) 


%% Set up GLM for each subject

% initialize_GLM_stuff(subjnum)
% 
% finalize_GLM_stuff(subjnum)

% Visualize GLM stuff for this subject

% anat = niftiRead('/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/CC/sess2/anat_T1_brain.nii');
% slice = 83;
% 
% GLM = niftiRead('/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/CC/GLMresults/stats.CC.nii.gz');
% stats1 = GLM.data(:,:,:,1,1);
% 
% figure
% imshow(anat.data(:,:,slice))
% hold on
% imagesc(stats1(:,:,slice,1,1)); colorbar


