function [ii_sess] = run_import_iEye(subj,subjs_with_wrong_freq)
%% run_import_iEye
% use iEye package to load edf, preprocess data, output usable data
% structure
ii_sess = []; %initialize as empty

% very first thing we want to do is define all parameters for processing
% (see ii_loadparams.m for default values)


addpath(genpath('../iEye-iEye_ts/')) %add outer path to grab iEye user functions
addpath(genpath('/System/Volumes/Data/d/DATA/hyper/software/eyelinkPoisonIvy/EYELINK/elcl/EDF/'))
setpref('iEye_ts','edf2asc_path','/Applications/EyeLink/') %EyeLink edf2asc pref for path
setpref('examples','edf2asc_path','/Applications/Eyelink')
%setpref('examples','edf2asc_path','/System/Volumes/Data/d/DATA/hyper/software/eyelinkPoisonIvy/EYELINK/elcl/EDF/') %EyeLink edf2asc pref for path
if sum(subjs_with_wrong_freq == subj)>0 %sometimes I forgot to set this correctly in the scanner room
    ifg_fn = 'p_1000hz.ifg';
    % account for that here
else
    ifg_fn = 'p_500hz.ifg'; %intended frequency of eye-tracking data
end


ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to
ii_params.valid_epochs =[1 2 3 4 5 6 7];
ii_params.trial_end_value = 7;   % XDAT value for trial end
ii_params.drift_epoch = [1 2 3]; % XDAT values for drift correction
ii_params.calibrate_epoch = 5;   % XDAT value for when we calibrate (feedback stim)
ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
ii_params.calibrate_mode = 'scale'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.plot_epoch = [2 3 4];  % what epochs do we plot for preprocessing?
ii_params.calibrate_limits = [2.5]; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)

ii_params.ppd = 33.6; % for scanner, 1280 x 1024 - convert pix to DVA

% all files we want to load are exfmri_r??.mat - so let's list them (we
% know they're in the same directory as this script)
tmp = mfilename('fullpath'); tmp2 = strfind(tmp,filesep);
root = tmp(1:(tmp2(end)-1));

files = dir(['data/subj' num2str(subj) '/eyetracking/']); 
filenames = string(char(files.name));
edf_files = filenames(contains(filenames,'.edf')&contains(filenames,'eowm')); 

edf_files = strrep(edf_files,' ',''); folder_name = string(char(files.folder));
temp = [edf_files(contains(edf_files,'sess1'),:); edf_files(contains(edf_files,'sess2'),:)];
edf_files = temp;

if ~isempty(edf_files) %some subjects have no eyetracking data
    edf_files = fullfile(folder_name(1:length(edf_files),:),edf_files);

    % create empty cell array of all our trial data for combination later on
    ii_trial = cell(length(edf_files),1);

    for ff = 1:length(edf_files)

        % what is the output filename?
        preproc_fn{ff} = strrep(edf_files{ff},'.edf','_preproc.mat');

        % run preprocessing!
        [ii_data, ii_cfg, ii_sacc] = ii_preproc(edf_files{ff},ifg_fn,preproc_fn{ff},ii_params);

        if ff == 1
            % plot some features of the data
            % (check out the docs for each of these; lots of options...)
            ii_plottimeseries(ii_data,ii_cfg); % pltos the full timeseries

            ii_plotalltrials(ii_data,ii_cfg); % plots each trial individually

            ii_plotalltrials2d(ii_data,ii_cfg); % plots all trials, in 2d, overlaid on one another w/ fixations
        end

        % score trials
        % default parameters should work fine - but see docs for other
        % arguments you can/should give when possible
        [ii_trial{ff},ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc);

    end

    ii_sess = ii_combineruns(ii_trial);
end


end