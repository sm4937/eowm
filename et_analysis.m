%% Run eyetracking analysis for EOWM
% what does this entail? I don't know. 
% I want to see increased pupil size in the hard versus the easy condition
% got new package that Grace seems to have worked on here:
%https://github.com/clayspacelab/iEye/tree/iEye_ts
function [t] = et_analysis(t,conditions,correct)
    subj = t.subj;
    files = dir(['data/subj' num2str(subj) '/eyetracking/']); 
    filenames = string(char(files.name));
    import_flag = sum(contains(filenames,'preproc_timeseries'))==0; %for now it's just whether there's any, should eventually do run #
    if import_flag
        ii_sess = run_import_iEye(subj); 
        save(['data/subj' num2str(subj) '/eyetracking/subj' num2str(subj) '_ii_sess.mat'],'ii_sess')
        close all
    else
        load(['data/subj' num2str(subj) '/eyetracking/subj' num2str(subj) '_ii_sess.mat'])
    end
    %function to run top level of example script in (/Users/sarah/Documents/MATLAB/iEye-iEye_ts/examples/example_analysis.m)
    % ii_sess has one cell per trial within each structure variable (including Pupil
    % size data)
    % so for 6 runs of 12 trials, there are 72 cells of information inside

    % for reference, XDAT TAGS:
    %1. Cue
    %2. Target
    %3. Delay
    %4. Test
    %5. Feedback
    %6. Post-feedback
    %7. ITI

    delaytag = 3; stimtag = 2; 
    window = 3000; %how many frames to analyze, of late or early data? 3000 this gives 6 seconds of samples
    fixation_mandatory = [1 2 3 4]; %how much deviation okay?  deg visual angle? (more than 5? how many DOV is the fixation?)
    %let's extract info from delay period (XDAT 3)
    
    total_excluded = 0;
    
    ntrials = length(ii_sess.Pupil);
    for tt = 1:ntrials %cycle over trials
        
        if ii_sess.n_sacc(tt) == 0 %no break from fixation
            % is this the right way to do this?
            relevant = ii_sess.XDAT{tt}==delaytag;
            delay_size = ii_sess.Pupil{tt}(relevant); %grab size of delay pupil
            contrast = nanmean(ii_sess.Pupil{tt-(tt>1)}(ii_sess.XDAT{tt-(tt>1)}==7)); %grab ITI period
             %compare to ITI pupil size from that trial before (unless trial 1)

            delay_size_avg(tt) = nanmean(delay_size)-contrast;
            delay_size_late(tt) = nanmean(delay_size(end-window:end))-contrast;
            delay_size_early(tt) = nanmean(delay_size(1:window))-contrast;
            delay_size_timecourse(tt,:) = NaN(1,6010); %initialize same size
            % instead of just subtracting contrast, compute a % signal change:
            % PSC_trial = 100 * ((delay_size)./contrast) - 1)
            % delay_size_timecourse(tt,1:length(delay_size)) = 100*((delay_size./contrast)-1);
            delay_size_timecourse(tt,1:length(delay_size)) = delay_size-contrast;

            relevant = ii_sess.XDAT{tt}==stimtag;
            stim_size = ii_sess.Pupil{tt}(relevant)'; %grab size of stimulus presentation pupil
            pres_size_avg(tt) = nanmean(stim_size)-contrast;
            pres_size_timecourse(tt,:) = NaN(1,260);
            % pres_size_timecourse(tt,1:length(stim_size)) = 100*((stim_size./contrast)-1);
            pres_size_timecourse(tt,1:length(stim_size)) = stim_size-contrast;

            no_breaks = sum(ii_sess.XDAT{tt}==fixation_mandatory,2)==1;
            
        else % they broke from fixation at some point (not being picky about this at the moment)
            
            % NaN it all out
            delay_size_avg(tt) = NaN;
            delay_size_late(tt) = NaN;
            delay_size_early(tt) = NaN;
            delay_size_timecourse(tt,:) = NaN(1,6010); %initialize same size

            pres_size_avg(tt) = NaN;
            pres_size_timecourse(tt,:) = NaN(1,260);
            
            total_excluded = total_excluded + 1;

        end % of saccade if statement
        
    end
    
    if ntrials > length(conditions) % if there's more pupil data than task data (uncommon)
        temp = NaN(ntrials,2); temp(1:length(conditions),:) = conditions(1:length(conditions),:); 
        temp_cor = NaN(ntrials,1); temp_cor(1:length(correct),:) = correct(1:length(correct),:);
        conditions = temp; correct = temp_cor;
    elseif length(conditions) > ntrials % if there's more task data than pupil data (common)
        conditions = conditions(1:ntrials,:); correct = correct(1:ntrials,:);
    end
    t.cond_pres_pupil_size = [nanmean(pres_size_avg((conditions(:,1))==1&(correct==1))) nanmean(pres_size_avg((conditions(:,1)==2)&(correct==1)))];
    t.cond_delay_pupil_size = [nanmean(delay_size_avg((conditions(:,1))==1&(correct==1))) nanmean(delay_size_avg((conditions(:,1)==2)&(correct==1)))];
    t.cond_late_pupil_size = [nanmean(delay_size_late((conditions(:,1)==1)&(correct==1))) nanmean(delay_size_late((conditions(:,1)==2)&(correct==1)))];
    t.cond_early_pupil_size = [nanmean(delay_size_early((conditions(:,1)==1)&(correct==1))) nanmean(delay_size_early((conditions(:,1)==2)&(correct==1)))];
    t.delay_tc_pupil_size = nanmean(delay_size_timecourse,1); %timecourse of delay period pupil size
    t.pres_tc_pupil_size = nanmean(pres_size_timecourse,1); %same as above for target dot presentation epoch
    t.hard_delay_tc_pupil_size_incorrect = nanmean(delay_size_timecourse((conditions(:,1)==2)&(correct==0),:),1);
    t.easy_delay_tc_pupil_size_incorrect = nanmean(delay_size_timecourse((conditions(:,1)==1)&(correct==0),:),1);
    t.hard_delay_tc_pupil_size = nanmean(delay_size_timecourse((conditions(:,1)==2)&(correct==1),:),1);
    t.easy_delay_tc_pupil_size = nanmean(delay_size_timecourse((conditions(:,1)==1)&(correct==1),:),1);  
    t.hard_pres_tc_pupil_size = nanmean(pres_size_timecourse((conditions(:,1)==2)&(correct==1),:),1);
    t.easy_pres_tc_pupil_size = nanmean(pres_size_timecourse((conditions(:,1)==1)&(correct==1),:),1);  %time course of pupil size for target dot presentation, correct trials only
    t.hard_pres_tc_pupil_size_incorrect = nanmean(pres_size_timecourse((conditions(:,1)==2)&(correct==0),:),1);
    t.easy_pres_tc_pupil_size_incorrect = nanmean(pres_size_timecourse((conditions(:,1)==1)&(correct==0),:),1);
    t.pct_excluded = total_excluded./ntrials;
    
end



