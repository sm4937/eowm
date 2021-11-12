%function [] = initialize_GLM_stuff(subjnum,session)
% initialize_GLM_stuff Get up the files necessary for running a GLM on BOLD
% data through AFNI (base code written by Tommy Sprague/Masih Rahmati)
clear all
subjnum = 4; sessions = 2;

load('subjinits.mat')
inits = subjinits{subjnum};

%% Do moving of files, renaming of some, into GLMresults folder on datb
% Parameters specified by do_befGLM.sh documentation written by Masih

% Set up ROOTDIR
datadir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults/data'];
mkdir(datadir);
stimdir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults/predictors'];
mkdir(stimdir);

for session = 1:sessions

    % Set the location of your raw data here
    raw_data_dir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/sess' num2str(session) '/'];
    preprocessed_files = dir(raw_data_dir); 
    
    if ~contains(ls(datadir),['surf_volreg_sess' num2str(session) '_detrend.r08'])
        % if this is true, this copying has already happened
        % skip all this laborious stuff

        raw_data_prefix = 'surf_volreg_detrend';
        motion_prefix = 'dfile';
        motion_folder = 'SEalign.results';
        raw_data_files = []; motion_folder_names = []; old_motion_files = []; old_motion_folders = [];

        for file = 1:length(preprocessed_files)
           name = preprocessed_files(file).name;
           if contains(name,raw_data_prefix)
               raw_data_files = [raw_data_files; name];
           end
           if contains(name,motion_folder)
               motion_folder_names = [motion_folder_names; name];
           end
        end

        for f = 1:size(raw_data_files,1)
            temp  = strrep(raw_data_files(f,:),'volreg_',['volreg_sess' num2str(session) '_']);
            % put session name in 
            runnum = num2str(f);
            if f < 10
                runnum = [num2str(0) runnum];
            end
            temp = strrep(temp,runnum,['.r' runnum]);
            new_data_files(f,:) = temp;
        end

        % Find names of exact files, relevant to running the GLM
        assignments = [];
        for f = 1:size(motion_folder_names,1)
            folder = motion_folder_names(f,:);
            %pat = '.(\w*).nii.gz';
            pat = [inits '_sess' num2str(session) '_r(\w*)to(\w*)_SEalign.results'];
            answer = regexp(folder, pat, 'tokens');
            assignments(f,1) = str2num(answer{1}{1}); assignments(f,2) = str2num(answer{1}{2}); 
            motion_folder_contents{f} = dir([raw_data_dir motion_folder_names(f,:) '/']);

            for file = 1:size(motion_folder_contents{f})
                name = motion_folder_contents{f}(file).name;
                if contains(name,motion_prefix)
                    try
                        old_motion_files = [old_motion_files; name];
                        old_motion_folders = [old_motion_folders; folder];
                    catch
                        % ignore .rall files for now
                    end
                end
            end

        end

        for file = 1:size(raw_data_files,1)
            copyfile([raw_data_dir raw_data_files(file,:)],[datadir '/' new_data_files(file,:)])
        end

        % Drag .1D motion files over for each run, and make sure they're properly
        % labeled
        for run = 1:size(raw_data_files,1)
            if run < 10
                new_motion_filename = ['dfile_sess' num2str(session) '.r0' num2str(run) '.1D'];
            else
                new_motion_filename = ['dfile_sess' num2str(session) '.r' num2str(run) '.1D'];
            end
            copyfile([raw_data_dir old_motion_folders(run,:) '/' old_motion_files(run,:)], ...
                [datadir '/' new_motion_filename]);
        end
%         new_motion_filename = ['dfile_sess' num2str(session) '.rall.1D'];
%         copyfile([raw_data_dir old_motion_folders(run,:) '/' old_motion_files(run,:)], ...
%                 [datadir new_motion_filename]);
        
    end % of skipping if statement
end % of looping over all sessions

%% Make text files for each predictor of interest

timing = get_precise_task_timing(subjnum,sessions);
predictors = [1 3 4]; 
% 1 : Cue & target stimulus presentation rolled into one
% 3 : Delay period
% 4 : Target/response period

% What are your six predictors?
% 1. cue/stim easy (t1_easy)
% 2. cue/stim hard (t1_hard)
% 3. delay easy (t3_easy)
% 4. delay hard (t3_hard)
% 5. response easy (t4_easy)
% 6. response hard (t4_hard)

nruns = max(timing.run(timing.sess==1)) + max(timing.run(timing.sess==2));
ntrials = max(timing.trial);
t_1_easy = -1*ones(nruns,ntrials+1); t_1_hard = t_1_easy; 
t_3_easy = -1*ones(nruns,ntrials+1); t_3_hard = t_3_easy;
t_4_easy = -1*ones(nruns,ntrials+1); t_4_hard = t_4_easy;

for sess = 1:max(timing.sess)
    for run = 1:max(timing.run)
        
        overallrun = run;
        if sess == 2
            overallrun = run + max(timing.run(timing.sess==1)); 
        end
        
        for trial = 1:ntrials
            
            this_trial = timing(timing.trial == trial & timing.run == run & timing.sess == sess,:);
            
            conds(overallrun,trial) = unique(this_trial.conds);
            
            for p = 1:length(predictors)
                cond = conds(overallrun,trial);
                if cond == 1
                    eval(['t_' num2str(predictors(p)) '_easy(overallrun,trial) = this_trial.times(this_trial.events==predictors(p));'])
                elseif cond == 2
                    eval(['t_' num2str(predictors(p)) '_hard(overallrun,trial) = this_trial.times(this_trial.events==predictors(p));'])
                end
            end
            
        end
    end
end

t_1_easy = round(t_1_easy,2); t_1_hard = round(t_1_hard,2);
t_3_easy = round(t_3_easy,2); t_3_hard = round(t_3_hard,2);
t_4_easy = round(t_4_easy,2); t_4_hard = round(t_4_hard,2);

writematrix(t_1_easy,[stimdir '/cuestim_easy.txt'])
writematrix(t_1_hard,[stimdir '/cuestim_hard.txt'])
writematrix(t_3_easy,[stimdir '/delay_easy.txt'])
writematrix(t_3_hard,[stimdir '/delay_hard.txt'])
writematrix(t_4_easy,[stimdir '/resp_easy.txt'])
writematrix(t_4_hard,[stimdir '/resp_hard.txt'])

rootdir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults'];

disp('GLM prep is finished!!!')
disp('Please run the following command in a terminal:')
disp(['cd ' rootdir])
disp('bash do_befGLM.sh \')
disp(['	-subject ' inits ' \'])
disp(['	-dataDIR ' datadir ' \'])
disp(['	-stimDIR ' stimdir ' \'])
disp('	 2>&1 | tee output.do_befGLM')

