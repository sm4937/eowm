function [t] = makeEOWMdatatable(subjnum)
%makeEOWMdatatable Load up each data file from the subject's folder,
% concatenate into one big table with everything in place
t = table;

folder_name = ['data/subj' num2str(subjnum)];
% addpath(folder_name)
files = dir(folder_name); files = char(files.name);
files = files(4:end,:); %get rid of . and .. at folder head
files(contains(string(files),'eyetracking'),:) = []; %remove eyetracking files from this stage
files(contains(string(files),'fmri'),:) = []; %also the fmri folder

for f = 1:size(files,1) %rows
    allp{f} = load([folder_name '/' files(f,:)]); %pull p's out of each file
end

%grab relevant information out of allps and put into t
rts = []; conditions = []; wm_ang_all = []; correct = [];
%initialize blank variables

for run = 1:length(allp)
    p = allp{run}.p; %load up everything from p
    correct = [correct; p.correct(:,p.run+1)];
    conditions = [conditions; p.conditions];
    deltas = p.deltas_all(1:end-1,:); % the final deltas_all is unused
    target_accuracy = p.target_accuracy;
    wm_ang_all = [wm_ang_all; p.wm_ang];
    if find(contains(fieldnames(p),'MGS_start'))>0
        p.test_start = p.MGS_start;
    end
    rts = [rts; p.feedback_start-p.test_start];
end

% use maxruns to set default length of run by run means, etc.
correct_vec = correct(:); rt_vec = rts(:);
cum_accuracy = cumsum(correct_vec)./(1:length(correct_vec))';
t.subj = subjnum;
t.overall_accuracy = cum_accuracy(end);
t.mean_rt = nanmean(rts(:));
t.cond_accuracy = [mean(correct_vec(conditions(:,1)==1)) mean(correct_vec(conditions(:,1)==2))];
t.cond_rt = [mean(rt_vec(conditions(:,1)==1)) mean(rt_vec(conditions(:,1)==2))];
cum_cond_accuracy(:,1:2) = [cumsum(correct_vec(conditions(:,1)==1))./(1:sum(conditions(:,1)==1))' ...
    cumsum(correct_vec(conditions(:,1)==2))./(1:sum(conditions(:,1)==2))'];
% grab accuracy by quadrant
quad_identity(wm_ang_all(:)<90) = 1;
quad_identity(wm_ang_all(:)<180&wm_ang_all(:)>90) = 2;
quad_identity(wm_ang_all(:)<270&wm_ang_all(:)>180) = 3;
quad_identity(wm_ang_all(:)<360&wm_ang_all(:)>270) = 4;
for q = 1:4
    t.accuracy_by_quad(:,q) = mean(correct(quad_identity==q));
    t.rt_by_quad(:,q) = nanmean(rts(quad_identity==q));
end
t.deltas = deltas(end,:);

% grab pupil information
subjs_with_wrong_freq = [10, 11];
% did I fix this issue correctly?
% TO DO: CHECK THIS ANALYSIS AGAIN
t = eyetracking_analysis(t,conditions,correct_vec,subjs_with_wrong_freq);

%check out the staircase
plot_flag = false;
if plot_flag
    figure
    subplot(1,2,1)
    plot(1:size(deltas,1),deltas(:,1),'b')
    hold on
    plot(1:size(deltas,1),deltas(:,2),'r')
    legend({'Easy','Hard','Easy target','Hard target'},'Location','Best')
    title(['Subj ' num2str(subjnum) ' deltas thru staircase'])
    
    subplot(1,2,2)
    plot(1:length(cum_accuracy),cum_accuracy,'k')
    hold on
    plot(1:length(cum_cond_accuracy),cum_cond_accuracy(:,1),'b')
    plot(1:length(cum_cond_accuracy),cum_cond_accuracy(:,2),'r')
    plot(xlim,[target_accuracy(1) target_accuracy(1)],'b--')
    plot(xlim,[target_accuracy(2) target_accuracy(2)],'r--')
    ylabel('p(correct)'); xlabel('trial #')
    title('Accuracy over trials')
    fig = gcf; fig.Color = 'w';
end


end
