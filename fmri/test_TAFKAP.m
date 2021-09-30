%% Test TAFKAP 
% Run a generate/recover on TAFKAP model fitting 
clear all
addpath(genpath('../../TAFKAP'))

%% Simulate data

nvoxels = 100;
ntrials = 120;
% 180 degrees of stimuli by 75 trials or whatever
stimuli = randperm(180,ntrials);
stim_width = 5;
screen = zeros(ntrials,180);
% Set up noisy stimulus presentation
for trial = 1:ntrials
    screen(trial,:) = normpdf(1:180,stimuli(trial),stim_width);
end
runs = reshape(repmat([1:10],12,1),ntrials,1);

RFs = zeros(50,180);
RF1 = normpdf(1:180,60,15);
RF2 = normpdf(1:180,120,15);
RFs = [repmat(RF1,50,1); repmat(RF2,50,1)];

BOLD = RFs*screen';
% One voxel per row, 1 trial per column
BOLD = BOLD';
% Now one trial per row, 1 voxel per column

test = false(ntrials,1);
test(runs==10) = true;
train = ~test;

p = struct;
p.stimval = stimuli;
p.runNs = runs;
p.test_trials = test;
p.train_trials = train;
p.stim_type = 'circular';

[est, unc, liks, hypers]  = TAFKAP_Decode(BOLD,p)

