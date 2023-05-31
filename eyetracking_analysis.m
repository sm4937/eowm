%% Run eyetracking analysis for EOWM
% THIS IS CURRENTLY A WORK IN PROGRESS
% THERE ARE ANALYSES WITH RAW PUPIL SIZE
% AND ANALYSES WITH PREPROCESSED/FILTERED PUPIL SIZE
% ALL MIXED IN TOGETHER
% IT'S A MESS!

% I hope to see increased pupil size in the hard versus the easy condition
% the new iEye package that this relies on is here:
%https://github.com/clayspacelab/iEye/tree/iEye_ts

function [t] = eyetracking_analysis(t,conditions,correct,subjs_with_wrong_freq)
subj = t.subj;
files = dir(['data/subj' num2str(subj) '/eyetracking/']);
filenames = string(char(files.name));

import_flag = sum(contains(filenames,'_ii_sess'))==0; %for now it's just whether there's any, should eventually do run #
if import_flag
    ii_sess = run_import_iEye(subj,subjs_with_wrong_freq);
    save(['data/subj' num2str(subj) '/eyetracking/subj' num2str(subj) '_ii_sess.mat'],'ii_sess')
    close all
else
    load(['data/subj' num2str(subj) '/eyetracking/subj' num2str(subj) '_ii_sess.mat'])
end
%function to run top level of example script in (/Users/sarah/Documents/MATLAB/iEye-iEye_ts/examples/example_analysis.m)
% ii_sess has one cell per trial within each structure variable (including Pupil
% size data)
% so for 6 runs of 12 trials, there are 72 cells of information inside

% extract saccade timing so we know when "illegal" saccades were made, as
% opposed to okay ones (during ITI)

ntrials = length(ii_sess.Pupil); 
saccade_flag = false(ntrials,1);

files = dir(['data/subj' num2str(subj) '/eyetracking/']);
filenames = string(char(files.name));
sacc_files = filenames(contains(filenames,'_preproc_sacc.mat'));
for f = 1:length(sacc_files)
    file = char(sacc_files(f,:));
    file = strrep(file,' ','');
    load(['data/subj' num2str(subj) '/eyetracking/' file])
    sess = regexp(file,'_sess(\d*)_','tokens');
    sess = str2num(sess{1}{1});
    run = regexp(file,'_run(\d*)_','tokens');
    run = str2num(run{1}{1})+1;
    
    epoch_start = ii_sacc.epoch_start;
    fixation_mandatory = [1 2 3 4];
    
    % flag saccades during cue, stim, delay, or test periods
    broken = sum(epoch_start==fixation_mandatory,2)>0 & ...
        ii_sacc.peakvelocity > 30 & ...
        ii_sacc.duration > 0.0075 & ...
        ii_sacc.amplitude > 0.25;
        
    %saccades classified in iEye with following criteria
    %ii_params.sacc_velocity_thresh = 30
    %ii_params.sacc_duration_thresh = 0.0075
    %ii_params.sacc_amplitude_thresh = 0.2500
        
    trial_start = ii_sacc.trial_start(broken);
    % get trial number based on where saccade criteria are met, and during
    % the wrong epoch
    real_trial_start = unique(trial_start) + (12*(run-1)) + (120*(sess-1));
    
    % filter out saccades which are during okay periods
    
    % get true trial number, accounting for n runs & n sessions
    saccade_flag(real_trial_start) = true;
    
end

% for reference, XDAT TAGS:
%1. Cue
%2. Target
%3. Delay
%4. Test
%5. Feedback
%6. Post-feedback
%7. ITI

sample_Hz = 500;
downsample_flag = false;
if sum(subj==subjs_with_wrong_freq)>0
    sample_Hz = 1000;
    downsample_flag = true; 
end

total_excluded = 0;
n_max_delay_samples = 6900;
n_max_stim_samples = 3500;
n_max_trial_samples = 8000;

% align the task and pupil data in their size
correct_unedited = correct;
if length(conditions) < ntrials
    addon = ntrials-length(conditions);
    conditions = [conditions; NaN(addon,2)];
    correct = [correct; NaN(addon,1)];
    correct_unedited = correct;
elseif length(conditions) > ntrials % if there's more task data than pupil data (common)
    conditions = conditions(1:ntrials,:); correct = correct(1:ntrials,:);
end

% initialize variables for trial-by-trial pupil size (delay period)
% analysis
delay_Y = []; delay_model = []; 
stim_Y = []; stim_model = [];
conditions_long_delay = []; trials_long_delay = [];
conditions_long_stim = []; trials_long_stim = [];
conditions_long_trial = []; trials_long_trial = [];
% I know the last variable has an unfortunate name
% It's going to contain a list of trial numbers which correspond to the
% full list of samples (bc there are many more samples than trials)
% And it's going to correspond to the full trial timecourse

% grab just delay, just stim, or whole timecourse
tags = [3 0 0 0 0; 2 0 0 0 0; 1 2 3 4 5];
labels = {'delay','stim','trial'};
all_delay_timecourses = [];
all_stim_timecourses = [];

for epoch = 1:2
    % cycle over delay and stimulus period,
    % always use full trial timecourse for baselining etc.
    
    tag = tags(epoch,:);
    Y = [];
    model = [];
    all_trial_timecourses = [];

    for tt = 1:ntrials %cycle over trials
        
        if ~saccade_flag(tt) %no break from fixation
            
            relevant = sum(ii_sess.XDAT{tt}==tags(3,:),2)>0;
            
            if tt == 1 % first trial has no preceding ITI period
                contrast = nanmean(ii_sess.Pupil{tt}(ii_sess.XDAT{tt}==1));
            else
                % use last "bit" - maybe last 1 second, or last 05 seconds -
                % of ITI instead of whole ITI - the variable
                % lengths are potentially obscuring some trial stuff
                ITI_pupil = ii_sess.Pupil{tt-1}(ii_sess.XDAT{tt-1}==7);
                ITI_pupil = ITI_pupil(end-sample_Hz:end);
                contrast = nanmean(ITI_pupil);
            end
            
            % get size of pupil for whole trial 
            pupil_size = ii_sess.Pupil{tt}(relevant) - contrast;
            %pupil_size = ii_sess.Pupil{tt}(relevant);
            % downsample this if the subject data was collected with
            % the wrong sampling frequency (double the frequency of
            % sampling)
            too_long_flag = (sum(relevant)>n_max_trial_samples);
            if downsample_flag & too_long_flag
                relevant(2:2:end) = false;
                pupil_size = ii_sess.Pupil{tt}(relevant);
            end
            
            % CLEAN DATA
            % (FILTER HIGH-FREQUENCY NOISE)
            % (INTERPOLATE NANs FROM BLINKS)
            size_filtered = clean_pupil_data(pupil_size,sample_Hz);
            % get this trial's data cleaned up
            all_trial_timecourses = [all_trial_timecourses; conditions(tt,1) transpose(size_filtered) NaN(1,abs(length(size_filtered)-n_max_trial_samples))];
            
            % % RETHINK WHEN YOUR CONTRAST COMES FROM!
            % ITI doesn't seem appropriate for delay period contrast
            % and during ITI period before the trial started
            if tt == 1 % first trial has no preceding ITI period
                contrast = nanmean(ii_sess.Pupil{tt}(ii_sess.XDAT{tt}==1));
            elseif epoch == 1 % delay period
                % use last "bit" - maybe last 1 second, or last 0.5 seconds -
                % of stim presentation period
                stim_pupil = ii_sess.Pupil{tt}(ii_sess.XDAT{tt}==2);
                stim_pupil = stim_pupil(end-100:end);
                contrast = nanmean(stim_pupil);
            end
            
            relevant = sum(ii_sess.XDAT{tt}==tag,2)>0;
            eval(['too_long_flag = (sum(relevant)>n_max_' labels{epoch} '_samples);'])
            if downsample_flag & too_long_flag
                relevant(2:2:end) = [];
            end
            
            filtered_within_epoch = size_filtered(relevant);
            Y = [Y; filtered_within_epoch];
            eval(['all_' labels{epoch} '_timecourses = [all_' labels{epoch} '_timecourses; ' ...
                'conditions(tt,1) transpose(filtered_within_epoch) ' ...
                'NaN(1,abs(length(filtered_within_epoch)-n_max_' labels{epoch} '_samples))];']);
            
            X_timecourse = interpnan_NT(ii_sess.X{tt}(relevant));
            Y_timecourse = interpnan_NT(ii_sess.Y{tt}(relevant));
            
            condition_code = [-0.5 0.5];
            try
                condition_for_model = condition_code(conditions(tt,1)).*ones(sum(relevant),1);
            catch
                condition_for_model = NaN(sum(relevant),1);
            end
            model = [model; ones(sum(relevant),1) condition_for_model contrast.*ones(sum(relevant),1) X_timecourse Y_timecourse];
            
            eval(['conditions_long_' labels{epoch} ' = [conditions_long_' labels{epoch} '; repmat(conditions(tt,1),sum(relevant),1)];'])
            eval(['trials_long_' labels{epoch} ' = [trials_long_' labels{epoch} '; repmat(tt,sum(relevant),1)];'])
            
        else % they broke from fixation at the wrong time, and in a substantive way
            
            total_excluded = total_excluded + 1;
            
        end % of saccade if statement
        
    end % of trial-by-trial loop
    
    % trim things that will mess with your regression (NaN's)
    Y(sum(isnan(model),2)>0,:) = [];
    eval(['conditions_long_' labels{epoch} '(sum(isnan(model),2)>0,:) = [];'])
    eval(['trials_long_' labels{epoch} '(sum(isnan(model),2)>0,:) = [];'])
    model(sum(isnan(model),2)>0,:) = [];
    
    eval([labels{epoch} '_Y = Y;'])
    %delay_Y = Y;
    
    %model is intercept, difficulty condition, baseline size, x position of
    %the eye, Y position of the eye
    
    % REGRESS OUT EYE MOVEMENTS, BASELINE SIZE
    eval(['[' labels{epoch} '_betas, ' labels{epoch} '_residuals] = deconvolve_pupil(Y,model);'])
    
end % of looping over epochs

%% FURTHER PREPROCESS PUPIL DATA

% PLOT WHAT THE CLEANED TIMECOURSE LOOKS LIKE NEXT TO THE ORIGINAL 
% OPTIONAL PLOT
plot_flag = false;
if plot_flag
    trial_length = sum(relevant);
    indices = 1:trial_length:length(delay_residuals);
    for i = 1:length(indices)-100
        plot(delay_residuals(indices(i):indices(i)+trial_length))
        hold on
        plot(zscore(delay_Y(indices(i):indices(i)+trial_length)),'k')
        hold off
    end
    clean_fig();
    title('Residual pupil size vs true data (one trial)')
end


valid_trials = unique(trials_long_delay);

% % want to model response to easy/hard conditions
t.delay_betas = delay_betas'; % pad with NaNs to standardize size across subjects

t.stim_betas = stim_betas';

trial_indices = false(1,length(correct_unedited));
trial_indices(valid_trials') = true; % no break from fixation
trial_indices = trial_indices & (correct_unedited==1)';
t.pupil_trial_indices = [trial_indices  false(1,240-length(trial_indices))];
% this variable can be used to index in to task variables (like condition
% vector) to align it with pupil data above (the beta values)

t.cleaned_pupilsize_by_trial{1} = all_trial_timecourses;
t.cleaned_pupilsize_all_delays{1} = all_delay_timecourses;
t.cleaned_pupilsize_all_stims{1} = all_stim_timecourses;
t.pct_excluded = total_excluded./ntrials;

end

%% code graveyard
% old way of running pupil GLM
% for ii = 1:length(valid_trials)
%     tt = valid_trials(ii);
%     cond = conditions(tt,1);
%     indices = ceil(linspace(1,sum(trials_long_delay==tt),7));
%     % break timecourses into 6 2-second chunks for stats analysis
% 
%     if correct(tt) == 1
%         
%         temp = NaN(1,n_max_delay_samples);
%         temp(1:sum(trials_long_delay==tt)) = delay_residuals(trials_long_delay==tt,:)';
%         delay_size_timecourse_all = [delay_size_timecourse_all; temp]; %overall matrix w/ all cleaned delay periods
%         
%         if cond == 1
% 
%             delay_size_timecourse_easy = [delay_size_timecourse_easy; temp];
%             binned_means_easy = [binned_means_easy; ...
%                 nanmean(temp(indices(1):indices(2))) nanmean(temp(indices(2):indices(3))) ...
%                 nanmean(temp(indices(3):indices(4))) nanmean(temp(indices(4):indices(5))) ...
%                 nanmean(temp(indices(5):indices(6))) nanmean(temp(indices(6):indices(7)))];
% 
%         elseif cond==2
% 
%             delay_size_timecourse_hard = [delay_size_timecourse_hard; temp];
%             binned_means_hard = [binned_means_hard; ...
%                 nanmean(temp(indices(1):indices(2))) nanmean(temp(indices(2):indices(3))) ...
%                 nanmean(temp(indices(3):indices(4))) nanmean(temp(indices(4):indices(5))) ...
%                 nanmean(temp(indices(5):indices(6))) nanmean(temp(indices(6):indices(7)))];
% 
%         end
%         
%         temp = NaN(1,n_max_stim_samples);
%         temp(1:sum(trials_long_stim==tt)) = stim_residuals(trials_long_stim==tt,:)';
%         stim_size_timecourse_all = [stim_size_timecourse_all; temp];
%         
%         temp = NaN(1,n_max_trial_samples);
%         temp(1:sum(trials_long_trial==tt)) = trial_residuals(trials_long_trial==tt,:)';
%         trial_timecourse_all = [trial_timecourse_all; temp];
% 
%     end % of correct if statement
%     
% end % of trial loop
% 
% ms_per_bin = 300;
% % determine bin length
% samples_per_bin = 500 * (1/1000) * ms_per_bin;
% % samples/second * second/millisecond = samples/milliseconds *
% % milliseconds = samples
% nbins = floor(length(trial_timecourse_all)./samples_per_bin);
% 
% %run sliding window average over all trials
% bin_means = [];
% for bin = 1:nbins
%     start_index = 1+(bin-1)*samples_per_bin;
%     bin_means = [bin_means (nanmean(trial_timecourse_all(:,start_index:(start_index+samples_per_bin)),2))];
% end
% 
% delay_betas = run_pupil_GLM(delay_size_timecourse_all,sample_Hz)
