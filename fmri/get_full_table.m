function [full] = get_full_table(subjnum,sessions)
% Align behavioral and fmri data in a big table

datapath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];
allbehavfiles = strsplit(ls(datapath),'.mat'); %get behavioral data to match up

TR = 0.75; % in seconds
long = [];

for sess = 1:sessions
    load([datapath '/fmri/sess' num2str(sess) '.mat']);
    % load up functional data for that session

    runs = length(funcdata);
    nTRs = size(funcdata{1}.data,4); %how many total TRs?

    for r = 1:runs
        run = r-1; long_run = [];
        behavidx = find(contains(allbehavfiles,['run' num2str(run) '_sess' num2str(sess)]));
        % find and load the behav file for this run
        beginidx = find(allbehavfiles{behavidx}=='o'); beginidx = beginidx(1)-1;
        behav = load([datapath '/' allbehavfiles{behavidx}(beginidx:end) '.mat']);
        times = 0:TR:(TR*nTRs); times = times(1:nTRs)';
        run_start = behav.p.trial_start(1);
        run_end = behav.p.iti_start(end)-run_start;

        all_times = [behav.p.trial_start behav.p.targ_start behav.p.delay_start behav.p.test_start behav.p.feedback_start behav.p.iti_start behav.p.iti_start+behav.p.itis];
        % make a column vector which starts at 0 and ends at the end of trial
        all_times = reshape([all_times-run_start]',numel(all_times),1);
        all_times = [0; behav.p.start_wait+all_times; all_times(end)+behav.p.end_wait];

        % mark all events that occur
        % 1 is trial start, 2 is target start, 3 is delay start
        % 4 is test start, 5 is feedback start, 6 is iti start, 7 is 
        % time betwen ITI end and next trial start
        all_events = [0; repmat([1; 2; 3; 4; 5; 6; 7],12,1); 0];
        all_conds = [0; reshape(repmat(behav.p.conditions(:,1),1,7)',84,1); 0];
        all_ang = [NaN; reshape(repmat(behav.p.wm_ang,1,7)',84,1); NaN];
        all_x_coords = [NaN; reshape(repmat(behav.p.wm_coords(:,1),1,7)',84,1); NaN];
        all_y_coords = [NaN; reshape(repmat(behav.p.wm_coords(:,2),1,7)',84,1); NaN];
        all_resps = [NaN; reshape(repmat(behav.p.resp,1,7)',84,1); NaN];
        all_trials = repmat([1:12],7,1); all_trials = [0; all_trials(:); 0];
        % 0 padding and extra seconds added on come from wait time at start
        % and end of each run

        for tt = 1:length(times)
            slice = times(tt);
            idx = find(all_times<slice);
            if ~isempty(idx)
                long_run = [long_run; subjnum sess r all_trials(idx(end)) slice all_events(idx(end)) all_conds(idx(end)) all_ang(idx(end)) all_x_coords(idx(end)) all_y_coords(idx(end)) all_resps(idx(end))];
            end
        end

        % get relative time in delay period to compare late delay to early
        % delay later on, etc.
        delay_times = zeros(size(long_run,1),1);
        for trial = 1:max(unique(long_run(:,4)))
            start = behav.p.start_wait + behav.p.delay_start(trial) - run_start;
            delay_times(long_run(:,6)==3&long_run(:,4)==trial,:) = long_run(long_run(:,6)==3&long_run(:,4)==trial,5) - start;
        end
        long_run = [long_run delay_times];
        long = [long; long_run];

    end
end % end of loop over session    
    
% label in table for clarity
full = table; 
full.subj = long(:,1);
full.sess = long(:,2);
full.run = long(:,3); 
full.trial = long(:,4);
full.time = long(:,5);
full.epoch = long(:,6); 
full.cond = long(:,7);
full.stimval = long(:,8);
full.x_target = long(:,9);
full.y_target = long(:,10);
full.resp = long(:,11);
full.delay_times = long(:,12);
full.overalltrial = cumsum([0; diff(full.trial)~=0]);
full.overallrun = cumsum([1; diff(full.run)~=0]);

% delay_times = zeros(height(full),1);
% full.delay_times = delay_times;

idxes = ([false; diff(full.overalltrial)~=0]);
eval(['task_subj' num2str(subjnum) ' = full(idxes,:);'])
save(['task_subj' num2str(subjnum) '.mat'], ['task_subj' num2str(subjnum)])

end

