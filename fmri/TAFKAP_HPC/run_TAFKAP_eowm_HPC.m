function [TAFKAP_output] = run_TAFKAP_eowm_HPC(data,task)
%run_TAFKAP_eowm Run Bayesian inverted encoding model (TAFKAP) through
%toolbox designed by 
addpath(genpath('TAFKAP'))

p = struct;
p.stim_type = 'circular';
%p.Nboot = 200;
ROI_list = fieldnames(data);
%ROI_list(contains(ROI_list,{'V2d','V2v','V3d','V3v','VO1','VO2','LO1','LO2','TO1','TO2'})) = [];

for rii = 1:length(ROI_list)
    disp(['starting at ROI ' num2str(rii) ': ' ROI_list{rii}])
    ROI_name = ROI_list(rii);
    eval(['BOLD = data.' ROI_name{1} ';'])
    % load the relevant BOLD data, ROI-specific, where
    % trial is row, and voxel is column
    
    if nansum(BOLD(:))~=0
    
        % THERE IS CURRENTLY AN ISSUE OF WHO HAS WHAT ROIS NOW

        % Make stimulus values go from 0 to 180
        temp = (task.stimval)./2; %BOLD((height(task)+1):end,:) = [];
        % Trim trials with no stimulus out
        toexclude = isnan(temp) & (sum(isnan(BOLD),2)>0);
        BOLD = BOLD(~toexclude,:); 
        p.stimval = temp(~toexclude);
        task_trimmed = task(~toexclude,:);
        p.runNs = task_trimmed.overallrun;

        % Trim voxels out with no BOLD (just NaNs)
        BOLD(:,sum(isnan(BOLD),1)>0) = [];

        %Z-score each voxel's time course, within each run
        for runii = 1:max(task_trimmed.overallrun)
            idx = task_trimmed.overallrun==runii;
            BOLD(idx,:) = zscore(BOLD(idx,:),0,1);
        end

        % Cross many-fold fitting, leave one run out
        temp_est = []; temp_unc = []; temp_liks = {}; temp_hyper = {};
        parfor fold = 1:max(task_trimmed.overallrun)
        %for fold = 1:max(task_trimmed.overallrun)
            thisp = p;
            thisp.test_trials = thisp.runNs == fold;
            thisp.train_trials = thisp.runNs ~= fold;
            % select the 1 run to leave out (to test)
            %trial_numbers = find(thisp.test_trials);

            % determine which rows are being tested, for later saving
            [this_est, this_unc, this_liks, this_hypers]  = TAFKAP_Decode(BOLD,thisp);

            %an ugly way to write this part, but suitable for parallel computing using parfor.
            %To improve, one can initialize these variables first. (Hsin wrote
            %this)
            temp_est(:,fold) = this_est;
            temp_unc(:,fold) = this_unc';
            %Pest(fold) = thisp;
            temp_liks{fold} = this_liks;
            temp_hyper{fold} = this_hypers;

        end %of folding & deocoding

        %Put decoded results into the shape we want
        est = []; unc = []; liks = []; hypers = [];
        for fold = 1:max(task_trimmed.overallrun)
            est(p.runNs == fold,:) = temp_est(:,fold);
            unc(p.runNs == fold,:) = temp_unc(:,fold);
            liks(p.runNs == fold,:) = temp_liks{fold};
            hypers(fold,:) = temp_hyper{fold};
        end

        eval(['TAFKAP_output.' ROI_name{1} '.est = est;'])
        eval(['TAFKAP_output.' ROI_name{1} '.unc = unc;'])
        % save results of decoding in subject-specific file
        eval(['TAFKAP_output.' ROI_name{1} '.liks = liks;'])
        eval(['TAFKAP_output.' ROI_name{1} '.hypers = hypers;'])
        % not sure how important this will be later on, but let's save it
        % all
        
    else % if BOLD data structure is empty
        
        eval(['TAFKAP_output.' ROI_name{1} '.est = NaN;'])
        eval(['TAFKAP_output.' ROI_name{1} '.unc = NaN;'])
        % save results of decoding in subject-specific file
        eval(['TAFKAP_output.' ROI_name{1} '.liks = NaN;'])
        eval(['TAFKAP_output.' ROI_name{1} '.hypers = NaN;'])
        % fill ROI struct with NaNs
        
    end
    
    save('TAFKAP_in_progress.mat','TAFKAP_output')

end % of looping over ROIs

delete('TAFKAP_in_progress.mat')
    
end % of function


