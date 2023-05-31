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
subjlist = [4 5 6 7 8 9 10 11]; 
nsessions = [2 2 2 2 2 2 2 2];
n = length(subjlist);

subjnum = 8;
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

save('current_full.mat','full','-v7.3')
% Note: To find PRF params for each subject, look for file names
% resembling:
% PRFparams = niftiread('/System/Volumes/Data/d/DATA/data/vRF_tcs/CC/RF1/CC_RF1_vista/RF_ss5_25mm-fFit.nii.gz');
% if you ever use y param from PRF params, apply -1 to it. it's a flipped
% output

%% Analyze functional data in terms of percent signal change,
% across all PRF-defined ROIs 


%load('current_full.mat')


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

for area = 2:size(hemi_list,1)
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
        
        if contains(ls(ROI_path_subj),ROI_filename)
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
            RF_thresh = 0.001;

            inRF = (normpdf(0:359,task.stimval,1) * RFs(VE>=selection_cutoff,:)')>RF_thresh;
            outRF = (normpdf(0:359,task.stimval,1) * anti_RFs(VE>=selection_cutoff,:)')>RF_thresh;
            % consider increasing this threshold
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
            % field voxelx
            tc_hard_inRF(sii,:) = nanmean(tc_hard_inRF_subj);
            tc_easy_inRF(sii,:) = nanmean(tc_easy_inRF_subj);
            tc_hard_outRF(sii,:) = nanmean(tc_hard_outRF_subj);
            tc_easy_outRF(sii,:) = nanmean(tc_easy_outRF_subj);

        end % end of cutting subjects out who don't have this ROI
        
    end % of cycling over participants

    %figure('Position',[0 0 1200 500])
    figure(1)
    subplot(3,2,area)
    % Plot time course of PSC for this ROI, differentiating between
    % conditions
    avg_tc_hard = nanmean(tc_hard(:,1:20)); avg_tc_easy = nanmean(tc_easy(:,1:20));
    SEM_hard = nanstd(tc_hard(:,1:20))/sqrt(n);
    SEM_easy = nanstd(tc_easy(:,1:20))/sqrt(n);

    errorbar(avg_tc_hard,SEM_hard,'Color',condcolors(1,:),'LineWidth',1.25,'DisplayName','Hard trials')
    xticks(0:2:37)
    xticklabels(0:TR*2:TR*37)
    hold on
    errorbar(avg_tc_easy,SEM_easy,'Color',condcolors(2,:),'LineWidth',1.25,'DisplayName','Easy trials')
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
    pause(5); saveas(fig,['PSC/Mean PSC.jpg']);
    
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
    pause(5); saveas(fig,'PSC/Condition differences PSC.jpg');
    %saveas(fig,['PSC/' figlabel ':' char(ROI_name{1}) '.jpg']);
    
    figure(3)
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
    pause(5); saveas(fig,['PSC/Mean PSC in and out RF timecourse.jpg']);
        
end % of cycling over areas of interest (ROIs)
