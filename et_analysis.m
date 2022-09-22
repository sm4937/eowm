%% Run eyetracking analysis for EOWM
% what does this entail? I don't know.
% I want to see increased pupil size in the hard versus the easy condition
% got new package that Grace seems to have worked on here:
%https://github.com/clayspacelab/iEye/tree/iEye_ts
function [t] = et_analysis(t,conditions,correct,subjs_with_wrong_freq)
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

% for reference, XDAT TAGS:
%1. Cue
%2. Target
%3. Delay
%4. Test
%5. Feedback
%6. Post-feedback
%7. ITI

sample_Hz = 500;

delaytag = 3; stimtag = 2;
window = 3000; %how many frames to analyze, of late or early data? 3000 this gives 6 seconds of samples
fixation_mandatory = [1 2 3 4]; %how much deviation okay?  deg visual angle? (more than 5? how many DOV is the fixation?)
%let's extract info from delay period (XDAT 3)

total_excluded = 0;
n_max_delay_samples = 6010;

ntrials = length(ii_sess.Pupil); 

Y = []; model = []; conditions_long = [];
% initialize variables for trial-by-trial pupil size (delay period)
% analysis

for tt = 1:ntrials %cycle over trials
    
    if ii_sess.n_sacc(tt) == 0 %no break from fixation
        
        % low-pass filtering SHOULD be applied on pupil data
        % ( to remove high frequency noise ) 
        relevant = ii_sess.XDAT{tt}==delaytag;
        delay_size = ii_sess.Pupil{tt}(relevant); %grab size of delay pupil
        relevant = ii_sess.XDAT{tt}==stimtag;
        stim_size = ii_sess.Pupil{tt}(relevant)'; %grab size of stimulus presentation pupil
        
        if length(delay_size)>n_max_delay_samples
            delay_size = delay_size(1:2:end);
            % downsample this if the subject data was collected with
            % the wrong sampling frequency (double the frequency of
            % sampling)
            stim_size = stim_size(1:2:end);
        end
        contrast = nanmean(ii_sess.Pupil{tt-(tt>1)}(ii_sess.XDAT{tt-(tt>1)}==7)); %grab ITI period
        %compare to ITI pupil size from that trial before (unless trial 1)
        
        % consider linearly interpolating pupil size instead of NaNing out
        % also consider visually inspecting all pupil/eye data for other
        % artifacts (not just blinks)
        
        delay_size_avg(tt) = nanmean(delay_size)-contrast;
        delay_size_late(tt) = nanmean(delay_size(end-window:end))-contrast;
        delay_size_early(tt) = nanmean(delay_size(1:window))-contrast;
        delay_size_timecourse(tt,:) = NaN(1,n_max_delay_samples); %initialize same size
        
        % instead of just subtracting contrast, compute a % signal change:
        % PSC_trial = 100 * ((delay_size)./contrast) - 1)
        delay_size_timecourse(tt,1:length(delay_size)) = 100*((delay_size./contrast)-1);
        %delay_size_timecourse(tt,1:length(delay_size)) = delay_size-contrast;
        
        
        pres_size_avg(tt) = nanmean(stim_size)-contrast;
        pres_size_timecourse(tt,:) = NaN(1,260);
        pres_size_timecourse(tt,1:length(stim_size)) = 100*((stim_size./contrast)-1);
        %pres_size_timecourse(tt,1:length(stim_size)) = stim_size-contrast;
        
        no_breaks = sum(ii_sess.XDAT{tt}==fixation_mandatory,2)==1;
        
        relevant = ii_sess.XDAT{tt}==delaytag;
        % TEMPORARILY REMOVING FIRST TWO SECONDS OF DATA PER TRIAL
        % BELIEVE TO BE STIMULUS-EVOKED PUPIL DILATION
        % WHEN WE RUN FULL TIMECOURSE ANALYSIS, THIS SHOULD BE REMOVED!!
        seconds_to_remove = 4;
        samples_to_remove = sample_Hz*seconds_to_remove;
        indices = find(relevant); relevant(indices(1:samples_to_remove)) = false;
        
        delay_size = ii_sess.Pupil{tt}(relevant); %grab size of delay pupil
        
        model = [model; ones(sum(relevant),1) contrast.*ones(sum(relevant),1) ii_sess.X{tt}(relevant) ii_sess.Y{tt}(relevant)];
        Y = [Y; delay_size];
        % new pupil timecourse 
        conditions_long = [conditions_long; repmat(conditions(tt),length(delay_size),1)];
        
    else % they broke from fixation at some point (not being picky about this at the moment)
        
        % NaN it all out
        delay_size_avg(tt) = NaN;
        delay_size_late(tt) = NaN;
        delay_size_early(tt) = NaN;
        % this is throwing an error because I think one of my subjects'
        % sessions was oversampled (i.e. sampled at 1000 Hz instead of
        % 500 Hz)
        delay_size_timecourse(tt,:) = NaN(1,6010); %initialize same size
        
        pres_size_avg(tt) = NaN;
        pres_size_timecourse(tt,:) = NaN(1,260);
        
        total_excluded = total_excluded + 1;
        
    end % of saccade if statement
        
end % of trial-by-trial loop



if ntrials > length(conditions) % if there's more pupil data than task data (uncommon)
    temp = NaN(ntrials,2); temp(1:length(conditions),:) = conditions(1:length(conditions),:);
    temp_cor = NaN(ntrials,1); temp_cor(1:length(correct),:) = correct(1:length(correct),:);
    conditions = temp; correct = temp_cor;
elseif length(conditions) > ntrials % if there's more task data than pupil data (common)
    conditions = conditions(1:ntrials,:); correct = correct(1:ntrials,:);
end

% run regression on pupil size on each timepoint of each trial (seconds 0 -
% 12 in the delay period)

delay_length_seconds = 12;
samples_per_second = 4;
% how many pupil timepoints you want per second, in this regression?
% basically, what level of granularity you interested in?
ntimepoints = delay_length_seconds*samples_per_second; %4 samples per second?
timepoints = linspace(1,n_max_delay_samples,ntimepoints);

%model = [ones(length(conditions(:,1)),1) conditions(:,1)==1 conditions(:,1)==2];

Y(sum(isnan(model),2)>0,:) = [];
conditions_long(sum(isnan(model),2)>0,:) = [];
model(sum(isnan(model),2)>0,:) = [];
% here we remove nans from the model & data (blinks)
% later on, we will linearly interpolate the data such that there are no
% NaNs to remove, and therefore more sound inferences can be made
betas = pinv(model)*Y;
% get effect of X,Y & baseline

Y_hat = model*betas;
% reconstruct pupil time course w/out X/Y or baseline information

residuals = Y - Y_hat;
% this is your new delay timecourse
residuals = zscore(residuals);

plot_flag = false;
if plot_flag
    trial_length = sum(relevant);
    indices = 1:trial_length:length(residuals);
    for i = 1:length(indices)-100
        plot(residuals(indices(i):indices(i)+trial_length))
        hold on
        plot(zscore(Y(indices(i):indices(i)+trial_length)),'k')
        hold off
    end
    clean_fig();
    title('Residual pupil size vs true data (one trial)')
end

t.cleaned_pupilsize_bycondition = [nanmean(residuals(conditions_long==1)) nanmean(residuals(conditions_long==2))];


% % want to model response to easy/hard conditions


% modelInv = pinv(model);
% for t = 1:length(timepoints)
%     time = timepoints(t);
% %    how far into the delay period are you?
%     data = delay_size_timecourse(:,time);
%     b(:,t) = modelInv * data;
% end
% 

% save stuff in usable, single-subject format (with common sizing to
% concatenate across subjects later)
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

t.alltrials_pres_mean_pupil_size = NaN(1,240);
t.alltrials_pres_mean_pupil_size(correct==1) = nanmean(pres_size_timecourse(correct==1,:),2);
% grab, on each trial, mean pupil size in entire stimulus presentation
% period

t.alltrials_delay_mean_pupil_size = NaN(1,240);
t.alltrials_delay_mean_pupil_size(correct==1) = delay_size_avg(correct==1);
t.alltrials_latedelay_mean_pupil_size = NaN(1,240);
t.alltrials_latedelay_mean_pupil_size(correct==1) = delay_size_late(correct==1);
t.alltrials_earlydelay_mean_pupil_size = NaN(1,240);
t.alltrials_earlydelay_mean_pupil_size(correct==1) = delay_size_early(correct==1);
% grab, on each trial, mean pupil size in late, early, entire delay
% period

end



