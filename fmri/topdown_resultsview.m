%% Preview group-level Top-Down analysis results
% Code written by Shanshan Li & Sarah Master May/June 2023

clear all

addpath(genpath('/Users/sarah/Documents/MATLAB/fmriTools'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/preproc_mFiles'))
addpath(genpath('/System/Volumes/Data/d/DATA/home/sarah/vistasoft_ts'))
addpath('../')
addpath(genpath('/Users/sarah/Documents/MATLAB/stats/'))
% pull up some pre-built functions

results = niftiRead('/System/Volumes/Data/d/DATB/datb/eowm_SM/TopDown_files/group_t&pV1_V2_V3.nii');
stats = results.data;
tmapn = stats(:,:,:,1);
pmapn = stats(:,:,:,2);

thresholdtmap = tmapn;
thresholdtmap(pmapn > 0.05) = NaN;
figure;
subplot(1,2,1)
vol3d('cdata', tmapn, 'texture', '3D');
view(3);
axis tight;
daspect([1 1 1]);
colormap jet;
colorbar;
title('tmap');
subplot(1,2,2)
vol3d('cdata', thresholdtmap, 'texture', '3D');
view(3);
axis tight;
daspect([1 1 1]);
colormap jet;
colorbar;
title('threshold, p < 0.05');

%% Visualize the data that produces the individual Z/group T values in Shanshan's analysis

condcolors = [190 0 110; 0 110 190]./255;

glmid = 'glm_TbT_del';

load('subjinits.mat')
subjnums = [4 5 6 7 8 10 11 12 13 14 16]; % subject number for EOWM
% subj 9 is left out here cause he doesn't have bilateral IPS2
N = length(subjnums);

%% Navigation map:
% 1. Pull beta matrix from rall_func_copy.nii.gz, which gives us trial by
% trial beta;

% 2. read task, which should have cond info on, figure out a way to read beta
% mat into two easy & hard maps

% 3. for different control regions, (and with V1_V2_V3 as seed region),
% assess visually & with mean Z value, the strength of the association

% cycle over 2 control regions of interest: sPCS and IPS2_3 (anterior
% parietal)

control_regions = {'sPCS','IPS2'};
for cii = 1:length(control_regions)
        
    figure(2+cii)
    for sii = 1:N
        
        subj = subjnums(sii);
        
        subjid = subjinits{subj};
        % define path, make new path
        glm_path = ['/d/DATC/datc/TopDown_SL/Mutha_GLM/eowm_data/' subjid '/GLMresults/' glmid '/'];
        correlation_path = [glm_path 'correlationResults/'];

        % Read TbT beta coefficients
        func_all_file = niftiRead([correlation_path 'rall_func_copy.nii.gz']);
        func_all_mat = func_all_file.data;
        
        %% Step 2: separate TbT_delay_coeff in to easy and hard
        % cond 1 = easy
        % cond 2 = hard
        
        % load stuff
        load(['task_subj' num2str(subj) '.mat'])
        eval(['task = task_subj' num2str(subj) ';'])
        task(ismember(task.cond, 0),:)=[]; %delete where cond = 0
        trials_easy = ismember(task.cond, 1);% find indeces
        trials_hard = ismember(task.cond, 2);
        % pick based on cond
        %     TbT_delay_coeff_easy = TbT_delay_coeff(:,:,:,trials_easy);
        %     TbT_delay_coeff_hard = TbT_delay_coeff(:,:,:,trials_hard);
        
        %% Step 3: separate correlation, with both acc and unc
        % load relevant stuff
        load(['TAFKAP_subj' num2str(subj) '_foveaROIs.mat'],'TAFKAP_output')
        
        % GET DECODING PERFORMANCE FROM SEED ROIs
        vis_decoding_acc = get_angular_distance(TAFKAP_output.V1_V2_V3.est.*2,task.stimval);
        
        vis_decoding_unc = TAFKAP_output.V1_V2_V3.unc;
                
        ROI_filename = ['/System/Volumes/Data/d/DATB/datb/eowm_SM/old_preprocessing_fmri/' subjinits{subj} '/ROIs/bilat.' control_regions{cii} '.nii.gz'];
        ROI = niftiRead(ROI_filename);
        ROI_delay_betas = niftiExtract(func_all_file,ROI);
        
        nvoxels = size(ROI_delay_betas,2);
        
        % Pick beta coefficients
        ntrials = (size(func_all_mat, 5) - (3 + 3 + 2))./3; %get number of trials
        % factoring in first 2 regressors are full F-stat, full R^2, then there are
        % 3 cue and 3 response regressors in the data
        delay_coeff_index = false(size(func_all_mat,5),1);
        delay_coeff_index(6:3:6+(3*(ntrials-1))) = true;
        
        % giant matrix of delay period coefficients
        TbT_delay_coeff = ROI_delay_betas(delay_coeff_index,:);
        % now TbT_delay_coeff is ntrials x N voxels

        
        subplot(4,4,sii+4)
        for vii = 1:nvoxels
            
            vector_easy = TbT_delay_coeff(trials_easy,vii);
            vector_hard = TbT_delay_coeff(trials_hard,vii);
            hold on
            scatter(vis_decoding_acc(trials_easy),vector_easy,[],condcolors(2,:),'Filled')
            scatter(vis_decoding_acc(trials_hard),vector_hard,[],condcolors(1,:),'Filled')
            
        end
        
        xlabel('Decoding error from V1-V2-V3')
        ylabel(['Delay period \beta from ' control_regions{cii}])
        title(['Subject ' num2str(sii)])
        xlim([0 180])

        clean_fig();
        hold off
        
        % correlation of betas and ROI results, trial-by-trial
        r_easy = NaN;
        r_hard = NaN;
        for vox = 1:nvoxels
            [r_easy(vox),p(vox)] = corr(TbT_delay_coeff(trials_easy,vox),vis_decoding_acc(trials_easy));
            [r_hard(vox),p(vox)] = corr(TbT_delay_coeff(trials_hard,vox),vis_decoding_acc(trials_hard));
        end
        
        %r-to-z transform
        zval_easy = atanh(r_easy);
        zval_hard = atanh(r_hard);
        
        subplot(4,1,1)
        hold on
        errorbar(sii+0.33,nanmean(zval_easy),nanstd(zval_easy)./sqrt(nvoxels),'LineWidth',2,'Color',condcolors(2,:),'DisplayName','Easy')
        errorbar(sii+0.66,nanmean(zval_hard),nanstd(zval_hard)./sqrt(nvoxels),'LineWidth',2,'Color',condcolors(1,:),'DisplayName','Hard')
        plot([sii+0.33 sii+0.66],[nanmean(zval_easy) nanmean(zval_hard)],'k--','LineWidth',2)
        
    end % of looping over subjects
    
    subplot(4,1,1)
    xticks((1:N)+0.5)
    ylabel('Mean Z value')
    xlabel('Subject')
    xticklabels(1:N)
    legend({'Easy','Hard'},'Location','Best')
    clean_fig();
    
end


