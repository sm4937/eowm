function [] = pull_fmri_data(subjnum,sess,datapath)

load('subjinits.mat','subjinits')
subject = subjinits{subjnum};

%datapath = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subject '/'];
addpath(genpath(datapath))
savepath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];

%anatpath = [datapath 'sub-' subject '/ses-anat/anat/sub-' subject '_ses-anat_desc-preproc_T1w.nii'];
%funcpath = [datapath 'sub-' subject '/ses-func/func/'];
anatpath = [datapath subject 'anat/'];
funcpath = [datapath 'sess' num2str(sess) '/'];


funcfiles = dir(funcpath); runs = [];
% From subject-specific folder, grab functional run numbers for later use
for f = 1:length(funcfiles)
    name = funcfiles(f).name;
    if contains(name,'surf_volreg_detrend')
        file = strsplit(name,'surf_volreg_detrend');
        runs = [runs; str2double(file{2}(1:2))];
    end
end

runs = unique(runs(~isnan(runs)));
for r = 1:length(runs)
    %functional_tag = ['sub-' subject '_ses-func_task-TASK_run-' num2str(r) '_space-T1w_desc-preproc_bold.nii'];
    if r < 10
        %functional_tag = ['surf_volreg_normPctDet0' num2str(r) '.nii.gz'];
        functional_tag = ['surf_volreg_detrend0' num2str(r) '.nii.gz'];  
    else
        functional_tag = ['surf_volreg_detrend' num2str(r) '.nii.gz'];
    end
    % This is an arbitrary functional file, probably not the one you
    % want to use, really
    filename = [funcpath functional_tag];
    [funcdata{r}] = niftiRead(filename);
end

%[anat] = niftiread(anatpath);
save([savepath '/fmri/sess' num2str(sess) '.mat'],'funcdata','-v7.3');
end

