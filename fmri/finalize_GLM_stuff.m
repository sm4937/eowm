%% Get GLM going, post-pre-processing of motion files & functional data
load('subjinits.mat')

subjnum = 4; sessions = 2;
inits = subjinits{subjnum};

% Set up ROOTDIR
datadir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults/data/'];
mkdir(datadir);
stimdir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults/predictors/'];
mkdir(stimdir);
destdir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults/data/procdata/'];
mkdir(destdir);
rootdir = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' inits '/GLMresults'];

    
filetag = 'scaled_surf_volreg_sess';
motiontag = 'dfile';

% Set the location of your raw data here
processed_files = dir(datadir); 

% if this is true, this copying has already happened
% skip all this laborious stuff

raw_data_files = []; motion_files = [];

for file = 1:length(processed_files)
   name = processed_files(file).name;
   if contains(name,filetag)
       raw_data_files = [raw_data_files; name];
   end
   if contains(name,motiontag)
       motion_files = [motion_files; name];
   end
end

for file = 1:size(raw_data_files,1)
    old_name = raw_data_files(file,:);
    pat = '_detrend.r(\w*).nii.gz';
    run = regexp(old_name, pat, 'tokens'); run = run{1}{1};
    pat = 'scaled_surf_volreg_sess(\w*)_detrend';
    sess = regexp(old_name, pat, 'tokens'); sess = str2num(sess{1}{1});
    overallrun = str2num(run) + (sess==2)*10;
    if overallrun < 10
        runtag = ['0' num2str(overallrun)];
    else
        runtag = num2str(overallrun);
    end
    new_name = strrep(old_name,['sess' num2str(sess) '_detrend.r' run],['detrend.r' runtag]);
    copyfile([datadir old_name],[destdir new_name]);
end

for file = 1:size(motion_files,1)
    old_name = motion_files(file,:);
    pat = '.r(\w*).1D';
    run = regexp(old_name, pat, 'tokens'); run = run{1}{1};
    pat = 'dfile_sess(\w*).r';
    sess = regexp(old_name, pat, 'tokens'); sess = str2num(sess{1}{1});
    overallrun = str2num(run) + (sess==2)*10;
    if overallrun < 10
        runtag = ['0' num2str(overallrun)];
    else
        runtag = num2str(overallrun);
    end
    new_name = strrep(old_name,['_sess' num2str(sess) '.r' run],['.r' runtag]);
    copyfile([datadir old_name],[destdir new_name]);
end

copyfile([datadir '/brainMask.nii.gz'],rootdir)

% All set! Short & sweet. Let's run the GLM:
disp('Run this command in your terminal:')
disp(['cd ' rootdir])
disp(['bash glm_6regressors.sh ' inits ' 2>&1 | tee output.glm.log'])
