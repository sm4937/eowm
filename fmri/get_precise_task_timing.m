function [timing_output] = get_precise_task_timing(subjnum,sessions)
% Align behavioral and fmri data in a big table

datapath = ['/Users/sarah/Documents/MATLAB/eowm/data/subj' num2str(subjnum)];

allbehavfiles = dir(datapath); filenames = [];
for file = 1:length(allbehavfiles)
    name = allbehavfiles(file).name;
    if contains(name,'eowm_v01')
        filenames = [filenames; name];
    end
end

timing_table = [];

for sess = 1:sessions
    session_files = filenames(contains(string(filenames),['sess' num2str(sess)]),:);
    runs = size(session_files,1);

    for r = 1:runs
        run = r-1; 
        behav = load([datapath '/' session_files(r,:)]);
        
        run_start = behav.p.trial_start(1);
        run_end = behav.p.iti_start(end)-run_start;

        all_times = [behav.p.trial_start behav.p.targ_start behav.p.delay_start behav.p.test_start behav.p.feedback_start behav.p.iti_start behav.p.iti_start+behav.p.itis];
        % make a column vector which starts at 0 and ends at the end of trial
        all_times = reshape([all_times-run_start]',numel(all_times),1);
        all_times = [0; behav.p.start_wait+all_times; all_times(end)+behav.p.end_wait];

        all_events = [0; repmat([1; 2; 3; 4; 5; 6; 7],12,1); 0];
        all_conds = [0; reshape(repmat(behav.p.conditions(:,1),1,7)',84,1); 0];
        all_trials = repmat([1:12],7,1); all_trials = [0; all_trials(:); 0];
        % 0 padding and extra seconds added on come from wait time at start
        % and end of each run
        
        timing_table = [timing_table; all_times all_events all_conds all_trials repmat(r,length(all_times),1) repmat(sess,length(all_times),1)];

    end
end % end of loop over session    

timing_output = table;
% label in table for clarity
timing_output.times = timing_table(:,1);
timing_output.events = timing_table(:,2);
timing_output.conds = timing_table(:,3);
timing_output.trial = timing_table(:,4);
timing_output.run = timing_table(:,5);
timing_output.sess = timing_table(:,6);

end

