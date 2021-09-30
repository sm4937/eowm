% eowm_task_v01.m
%
% 8/29/2018 - TCS & GH, adapted from MR code
% 01/05/2020 - SM - co-opted code structure from distractor task to make
% eowm version 1
%
% 
% TIMING INFO:
%
%
% XDAT INFO:
% Xdat tags like eye tracking triggers.
% Tags are as follows:
% 1. Cue
% 2. Target
% 3. Delay
% 4. Test
% 5. Feedback
% 6. Post-feedback
% 7. ITI

function eowm_v01(subj,run,scanner,session)

%try
p.expt_name = 'eowm_v01';

p.do_et = 1;
p.session = session;

p.TR = 3; % 4x multiband, so measured TR is 0.75, but "TR" for stim is 3

p.subj = subj;
p.run = run;
p.total_runs = 10;
p.scanner = scanner;
if ~exist(['./data/subj' num2str(subj)],'dir')
    mkdir(['./data/subj' num2str(subj)]);
    mkdir(['./data/subj' num2str(subj) '/eyetracking']);
end

p.filename = ['data/subj' num2str(p.subj) '/' p.expt_name '_subj' num2str(p.subj) '_run' num2str(p.run) '_sess' num2str(p.session) '_date' datestr(now,30)];
if p.do_et == 1
    p.eyedatafile = ['s' num2str(subj) 'r' num2str(run)];
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);

% define geometry of WM stimuli
% ------ size of relevant stim features, etc ------ %
p.wm_ecc = 10;     % deg [in behavioral room, max ecc of circle 15 deg]
p.cue_size = 0.50; % deg, was 0.55 changed by sarah for debugging
p.wm_size = 0.59;  % deg, size of WM dots
% was 0.65, changed by Sarah for debugging

% how far to be away from meridians (target pos); polar angle deg
p.merid_sep = 5; % must be this far (polar angle deg) away from vert/horiz

p.aperture_size = 15; % [or max ecc?]

%p.n_pos = 16; % and we'll use 2 offset angles

p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;

p.fix_size_mult = 1.25; % for pre/post experiment, fix is bigger
% new variable added by sarah, the *2 was throwing all sorts of dot sizing
% errors on home computer
p.fix_size_constant = 1.4; %used to be 2

% ------ color of relevant stim features ----------- %
% TODO: adjust??
p.bg_color  = 60;%20*[1 1 1];
p.fix_color = 200;%%75*[1 1 1];%[150 150 150];        % during trial/delay/etc

p.wm_color = p.fix_color;

p.go_color = p.fix_color;% 130*[1 1 1];%[255 255 255]; % when subj should choose, color of fix

p.dim_amt = 0.9; % multiply fix_color by this during ITI, start, end periods

% magenta:cyan RGB ratio: 2.46/1.18 (cyan must be 1.18/2.46x as intense as
% magenta)
p.cue_colors = [1 0 1;0 1 1]; % TODO: equiluminant! (magenta = no dist, cyan = dist)
p.cue_rel_lum = 0.5; % relative intensity, both colors scaled by this; magenta is this*255, cyan is this*1.18/2.46*255
p.cue_colors(1,:) = round(p.cue_colors(1,:)*p.cue_rel_lum*255);
p.cue_colors(2,:) = round(255*p.cue_rel_lum*(1.18/2.46)*p.cue_colors(2,:));
p.cue_colors = p.cue_colors(randperm(2,2),:);

p.feedback_colors = round([0 0.44*255 0;255 0 0;.305*255 .305*255 0]); % green, red, yellow for corr, incorr, miss

% define geometry, dot properties, etc of distractors

p.dist_jitter = 12; % ABSOLUTE randomize distractor position by + or - this (polar ang deg)

% ------ conditions ------ %

p.p_hard = 0.5; % probability of hard trial
% generate blank p.conditions that's this big, then consider blanks
% no-distractor
p.ntrials = 12; % # of trials per run

% column 1: easy/hard (1/2)
% column 2: counter/clockwise test (1/2)
p.conditions = [ones(p.ntrials*p.p_hard*0.5,1) ones(p.ntrials*p.p_hard*0.5,1); ...
    ones(p.ntrials*p.p_hard*0.5,1) 2*ones(p.ntrials*p.p_hard*0.5,1)];
p.conditions = [p.conditions; p.conditions(:,1)*2 p.conditions(:,2)];
% counterbalance probability of left/right so that it's always 50-50 in
% both easy and hard trials
p.target_accuracy = [0.9, 0.7];
p.correct_points = 1; p.incorrect_points = 0;

% shuffle these
p.rnd_idx = randperm(p.ntrials);
p.conditions = p.conditions(p.rnd_idx,:);
p.deltas_all = [30, 7]; %distance from target to test in easy, hard conditions (first trial)

% generate list of WM positions (all around the circle, 360, but no vertical or horizontal meridians)
% qtmp = floor(4*rand(p.ntrials,1)); % first, pick a quadrant (any quadrant)
% atmp = p.merid_sep + (90-2*p.merid_sep)*rand(p.ntrials,1); % randomize pos w/in quadrant
% % % above in polar ang
% p.wm_ang = qtmp*90+atmp; clear qtmp atmp;
possible_ang = 1:359;
possible_ang(mod(possible_ang,90)==0) = []; %no meridians
p.wm_ang_all = possible_ang(randperm(length(possible_ang),p.ntrials*p.total_runs*2)); %shuffle order
p.wm_ang_all = reshape(p.wm_ang_all,p.ntrials,p.total_runs*2); %create enough for 2 sessions per subject

%%
%p.TR = 3;
 
% ------ timing of trial events --------- %
% ITIs should be int_num * TR + 0.9 + 0.9 so that beginning of delay is
% locked to TR
 
% THINGS THAT HAPPEN BEFORE DELAY
p.cue_dur  = 1.0; % pre-cue: whether it's a hard or easy trial
p.targ_dur = 0.5; % WM target (every trial)

% Standardize delay to be always 12 seconds
p.delay_dur = 12;
 
% THINGS THAT HAPPEN AFTER DELAY (end of delay also locked to TR)
p.resp_dur = 0.8; % response period (time to make left/right response)
p.feedback_dur = 1.0; % feedback stimulus presentation- location of original stimulus in red/green 
p.post_feedback_dur = 1.0; %in seconds, the ITI
% try to be in range of ~ 6-12 s (avg 9)

% short/medium/long ITIs...
ITI_TRs = [1*ones(floor(p.ntrials/4),1); 2*ones(2*ceil(p.ntrials/4),1); 3*ones(floor(p.ntrials/4),1)];

p.itis = p.TR-(p.targ_dur+p.cue_dur) + ...
         p.TR*ITI_TRs  + ...
     mod(p.TR-(p.resp_dur+p.feedback_dur+p.post_feedback_dur),p.TR);

p.itis = p.itis(randperm(length(p.itis)));     

% -------- timing of experiment events ------- %
p.start_wait = 1 * p.TR + p.TR-(p.cue_dur+p.targ_dur);%2.4 + 0.9; % after first trigger [after dummys]
p.end_wait = 1*p.TR+mod(p.TR-(p.resp_dur+p.feedback_dur+p.post_feedback_dur+p.itis(1)),p.TR);

p.trial_dur = p.cue_dur + p.targ_dur + p.delay_dur + ...
              p.resp_dur + p.feedback_dur + p.post_feedback_dur + p.itis;

p.exp_dur = p.start_wait + p.end_wait + sum(p.trial_dur);

% being a bit redundant here...
p.trial_start = nan(p.ntrials,1);
p.targ_start  = nan(p.ntrials,1);
p.delay_start = nan(p.ntrials,1);
p.test_start   = nan(p.ntrials,1);
p.feedback_start = nan(p.ntrials,1);
p.iti_start   = nan(p.ntrials,1);
p.trial_end   = nan(p.ntrials,1);
p.behind_by = nan(p.ntrials,1); % keep track of timing errors, they should be basically tiny


% ------- keyboard stuff --------------------------- %
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    try
        p.esc_key = KbName('esc');
    catch
        p.esc_key = KbName('escape');
    end
end
p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.space = KbName('space');
p.resp_keys = [KbName('1!') KbName('2@')]; % LEFT, RIGHT

% ------- Screen setup, optics --------- %

% TODO: be a bit more clever here so that we can use this for practice
% outside scanner too....
if p.scanner == 0
    Screen('Preference', 'SkipSyncTests', 1) %changed from 21 on 3/12 
    p.desired_resolution = [1920 1080]; % desired resolution, compare to the 'actual' resolution below [this is for mbp laptop]
    p.desired_refresh_rate = 60;  % change this to 120 313
    % set to 60 for home computer testing by sarah, originally was 120
else
    Screen('Preference', 'SkipSyncTests', 1)
    p.desired_resolution = [1280 1024];% [1920 1080]; % % scanner
    p.desired_refresh_rate = 120; %changed from 120  on 312
end

FlushEvents('KbName')
% TODO: assert correct refresh rate, resolution; bail otherwise

tmp_res = Screen('Resolution',max(Screen('Screens')));
p.resolution = [tmp_res.width tmp_res.height];
p.refresh_rate = tmp_res.hz;
p.screen_height = 36.3; % cm

if sum(p.resolution==p.desired_resolution)~=2 || p.refresh_rate~=p.desired_refresh_rate
    sprintf('Unexpected resolution/RR: expected %i, %i, %i; found %i, %i, %i',p.desired_resolution(1),p.desired_resolution(2),...
        p.desired_refresh_rate,p.resolution(1),p.resolution(2),p.refresh_rate)
    error('spDist_scanner:displayError','Unexpected resolution/RR: expected %i, %i, %i; found %i, %i, %i',p.desired_resolution(1),p.desired_resolution(2),...
        p.desired_refresh_rate,p.resolution(1),p.resolution(2),p.refresh_rate);
end

p.screen_width = p.screen_height * p.resolution(1)/p.resolution(2); %
p.viewing_distance = 63; % cm
p.screen_height_deg = 2*atan2d(p.screen_height/2,p.viewing_distance);
p.screen_width_deg  = 2*atan2d(p.screen_width/2, p.viewing_distance);
p.ppd = p.resolution(2)/p.screen_height_deg;  % used to convert rects, positions later on

p.center = p.resolution/2;  % could do offset centers, etc?
% % PTB SCREEN OPENED HERE % %
%[windowPtr,rect]=Screen(‘OpenWindow’,windowPtrOrScreenNumber [,color] [,rect][,pixelSize][,numberOfBuffers][,stereomode][,multisample][,imagingmode][,specialFlags][,clientRect][,fbOverrideRect][,vrrParams=[]]);[windowPtr,rect]=Screen(‘OpenWindow’,windowPtrOrScreenNumber [,color] [,rect][,pixelSize][,numberOfBuffers][,stereomode][,multisample][,imagingmode][,specialFlags][,clientRect][,fbOverrideRect][,vrrParams=[]]);
if subj ~= 99 %not debugging
    [w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0]); HideCursor;
elseif subj == 99 %want small screen for debugging
    [w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0],[0 0 600 500]);
end
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);

%ListenChar(2); %turn matlab keyboard input off
p.ifi = Screen('GetFlipInterval',w);

% ------------------ generate rects we use later ---------------------- %
p.aperture_rect = CenterRectOnPoint([0 0 2 2] * p.ppd * p.aperture_size,p.center(1),p.center(2));
p.fix_rect_out  = CenterRectOnPoint([0 0 2 2] * p.ppd * p.fix_size_out, p.center(1),p.center(2));
p.fix_rect_in   = CenterRectOnPoint([0 0 2 2] * p.ppd * p.fix_size_in,  p.center(1),p.center(2));
p.aperture_radius = abs(p.aperture_rect(1)-p.aperture_rect(3))/2;%sqrt(sum(p.fix_rect_in(1:2)-p.aperture_rect(1:2)).^2); %inner x-value minus outer x-value
% --------- eyetracking ----------- %
if p.do_et == 1
    
    if p.scanner == 1
        Eyelink('SetAddress','192.168.1.5')
    end
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=200;
    el.calibrationtargetsize=2*p.wm_size;
    el.calibrationtargetwidth=1;
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
     % 192.168.1.5
    % sca
    EyelinkUpdateDefaults(el);

    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
   % SCANNER: right eye!!!!!!
    Eyelink('command','calibration_type=HV13'); % updating number of callibration dots
    s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate = 500');
    s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    
    
    % make sure that we get gaze data from the Eyelink
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_rect);
if run == 0 %if practice run, display instructions first
    runInstructions(w,p);
    p.all_conditions = zeros(0);
else %practice run already ran
    % load up old information about correctness, deltas, conditions
    filename = [p.expt_name '_subj' num2str(p.subj) '_run' num2str(p.run-1)];
    folder = ls(['data/subj' num2str(p.subj)]);
    specific_filename = folder(contains(string(folder),filename),:);
    old_p = load(['data/subj' num2str(subj) '/' specific_filename]); old_p = old_p.p;
    p.correct = old_p.correct; p.all_conditions = old_p.all_conditions; %keep list of previous behavior, etc.
    p.deltas_all = old_p.deltas_all; p.cue_colors = old_p.cue_colors; %keep counterbalancing of hard/easy colors/deltas
    p.wm_ang_all = old_p.wm_ang_all; %retain shuffle of WM target positions
end
%p.wm_ang = p.wm_ang_all(:,p.run+1+10*(p.session==2)); %grab run-specific target locations, adjusting for session #
p.wm_ang = p.wm_ang_all(:,p.run+1);
p.wm_coords = p.wm_ecc .* [cosd(p.wm_ang) sind(p.wm_ang)];

% draw up a message about which color is which difficulty level
cond_colors = p.cue_colors(1,1)>0; %magenta is first color, cyan second color
if cond_colors
    cond_color_names = {'Pink','Blue'};
else
    cond_color_names = {'Blue','Pink'};
end
this_message = ['Remember: ' cond_color_names{2} ' - more precise information.'];
this_message_2 = [cond_color_names{1} ' - less precise information.'];

Screen('FillRect',w,[0 0 0]);
draw_aperture();
draw_fixation(p,w);
DrawFormattedText(w,this_message,'center',p.center(2)-2*p.ppd,p.cue_colors(2,:));
DrawFormattedText(w,this_message_2,'center',p.center(2)+2*p.ppd,p.cue_colors(1,:));
Screen('Flip',w);

% check for esc, space.... 
resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        if p.do_et == 1
            Eyelink('ShutDown');
        end
        return;
    end
end
clear resp;
p.expt_start = GetSecs;

% blank screen - fix back to normal size
Screen('FillRect',w,[0 0 0]);
draw_aperture(); 
draw_fixation(p,w);
Screen('Flip',w);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 49);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


% ------ initial wait time (check for esc) ---------

resp = 0;
while (GetSecs-p.expt_start) < p.start_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressed during post-trigger wait time\n');
        ShowCursor;
        if p.do_et == 1
            Eyelink('ShutDown');
        end
        return;
    end
end
clear resp;

for tt = 1:p.ntrials
    % set up stimuli to be displayed HERE %
    % so, conditions = 2 means the test dot needs to be clockwise of the
    % 1st dot (negative change in angle)
    % counter-clockwise (1) means the test dot needs to be at a larger
    % angle
    if p.conditions(tt,2) == 2 %there's maybe a smarter way to do this but I don't want to think about it right now
        sign = -1;
    elseif p.conditions(tt,2) == 1
        sign = 1; %positive
    end
    delta = p.deltas_all(end,p.conditions(tt,1)); 
    %sqrt(x^2 + y^2) = delta ... be more sophisticated than this?
    p.test_ang(tt,:) = p.wm_ang(tt,:)+(delta*sign); %move - for clockwise, + for counterclockwise
    p.test_coords(tt,:) = p.wm_ecc .* [cosd(p.test_ang(tt,:)) sind(p.test_ang(tt,:))];
    
    % trial type cue (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % NOTE: I'm making this a "trial should have started" marker, to keep
    % timing from accumulating errors & to keep below code more readable
    trial_start = p.expt_start + p.start_wait + sum(p.trial_dur(1:(tt-1)));
    p.trial_start(tt) = GetSecs;
    p.behind_by(tt) = p.trial_start(tt)-trial_start;
    
    if p.do_et == 1

        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        Eyelink('Message','xDAT %i',1);
        
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        
    end

    tic
    while GetSecs < trial_start + p.cue_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        % fixation
        draw_fixation(p,w); 
        % cue (filled fix)
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*p.fix_size_constant-p.fix_pen,p.cue_colors(p.conditions(tt,1),:),p.center,2);         
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*p.fix_size_constant,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    
    % targets (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    p.targ_start(tt) = GetSecs;
    
    
    if p.do_et == 1
        Eyelink('Message','TargX %s', num2str(p.wm_coords(tt,1)));
        Eyelink('Message','TargY %s', num2str(p.wm_coords(tt,2)));
        Eyelink('Message','xDAT %i',2);
    end
    
    % show target
    while GetSecs < trial_start + p.cue_dur + p.targ_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        
        % fixation
        draw_fixation(p,w);
        Screen('Flip',w);
        
        if isnan(p.targ_start(tt))
            p.targ_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % Delay 1 (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);    
    end
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % fixation
        draw_fixation(p,w);
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    clear resp_start_time;
    
    % Response, display second stimulus(XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','TestX %s', num2str(p.test_coords(tt,1)));
        Eyelink('Message','TestY %s', num2str(p.test_coords(tt,2)));
        Eyelink('Message','xDAT %i',4);
    end

    while (GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay_dur + p.resp_dur) & resp == 0
        % wait until time expires OR response made, move on after response
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        draw_fixation(p,w);
        % draw second cue
        Screen('DrawDots',w,p.ppd*[1;-1].*p.test_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
    
        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*p.fix_size_constant,p.bg_color,p.center,2); 
        
        Screen('Flip',w);
        
        if isnan(p.test_start(tt))
            p.test_start(tt) = GetSecs;
        end
        
        % check for esc or response key
        [resp] = checkForResp(p.resp_keys, p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        p.resp(tt) = sum(find(p.resp_keys==resp)); %0 for non-response, 1 for 1 key, etc.
    end
    p.resp(tt)
    p.conditions(tt,2)
    % feedback (XDAT 5, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','TargX %s', num2str(p.wm_coords(tt,1)));
        Eyelink('Message','TargY %s', num2str(p.wm_coords(tt,2)));
        Eyelink('Message','xDAT %i',5); 
    end
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay_dur + p.resp_dur + p.feedback_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        % fixation
        draw_fixation(p,w);
        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*p.fix_size_constant,p.bg_color,p.center,2); 
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        if p.resp(tt)==0
            % draw yellow (feedback_colors 3)
            this_color = p.feedback_colors(3,:);
            this_message = 'Respond faster!';
        else
            if p.resp(tt) == p.conditions(tt,2)
                % draw green (1)
                this_color = p.feedback_colors(1,:);
                this_message = ['Correct! +' num2str(p.correct_points)];
            else
                this_color = p.feedback_colors(2,:);
                this_message = ['Incorrect! +' num2str(p.incorrect_points)];
            end
        end
        %Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*p.fix_size_constant-p.fix_pen,this_color,p.center,2);        
        % target dot, in feedback color
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, this_color, p.center, 2);
        % test dot
        Screen('DrawDots',w,p.ppd*[1;-1].*p.test_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        DrawFormattedText(w,this_message,'center',p.center(2)-1*p.ppd,this_color);
        Screen('Flip',w);
        
        if isnan(p.feedback_start(tt))
            p.feedback_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
 
    end
    
    % Post-feedback wait (XDAT 6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',6);    
    end
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay_dur + p.resp_dur + p.feedback_dur + p.post_feedback_dur 
    % don't change the screen just keep it the same
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % ITI (XDAT 7) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',7);
    end

    % run ITI with blank screen
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay_dur + p.resp_dur + p.feedback_dur + p.post_feedback_dur + p.itis(tt)
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        draw_fixation(p,w);        
        Screen('Flip',w);
        
        if isnan(p.iti_start(tt))
            p.iti_start(tt) = GetSecs;
            % save [note: in scanner, do this at beginning of ITI after first flip]
            save(p.filename,'p');
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor; ListenChar(1);

            if p.do_et==1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    trial_time = toc %make sure trial time still standardized
    
    % catalog trial ending time
    p.trial_end(tt) = GetSecs;
    % figure out staircase for next trial
    % DO STAIRCASING!
    p.correct(tt,run+1) = p.resp(tt)==p.conditions(tt,2); %did they press 1 or 2?
    %p.n_misses = sum(isnan(p.resp(p.conditions(:,1)==2)));
    p.acc = nanmean(p.correct(:));
    p.accuracy_run(run+1) = nanmean(p.correct(:,run+1));
    p.nmisses = sum(p.resp==0);
    % grab amount of distance between target and second cue
    p.all_conditions(end+1,:) = p.conditions(tt,:);
    relevant_correct = p.correct(:); relevant_correct = relevant_correct(p.all_conditions(:,1)==p.conditions(tt,1));
    new_deltas = p.deltas_all(end,:);
    new_deltas(p.conditions(tt,1)) = run_staircase(relevant_correct,p.target_accuracy(p.conditions(tt,1)),p.deltas_all(end,p.conditions(tt,1)));
    p.deltas_all(end+1,:) = new_deltas;
     %staircase result
    
end % end of trial loop

% ------- wait for p.end_wait -------- %
end_tmp = GetSecs;

resp = 0;
while (GetSecs-end_tmp) < p.end_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-experiment wait time\n');
        ShowCursor;ListenChar(1);
        return;
    end
end
clear resp;
p.end_expt = GetSecs;
p.points_earned = sum((p.correct(:)==1)*p.correct_points + (p.correct(:)==0)*p.incorrect_points);
save(p.filename,'p');

% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',50);
end

Screen('FillRect',w,[0 0 0]);
draw_aperture();
txt = sprintf('End of run %i',p.run);
txt2 = sprintf('Accuracy: %0.02f%%, %i Points Earned!',p.accuracy_run(run+1)*100,p.points_earned);
DrawFormattedText(w,txt, 'center',p.center(2)-5*p.ppd,p.fix_color);
DrawFormattedText(w,txt2,'center',p.center(2)-2*p.ppd,p.fix_color);
%fixation
draw_fixation(p,w); 
Screen('Flip',w);

if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    p.eyedatafile_renamed = [strrep(p.filename,['data/subj' num2str(subj) '/'],['data/subj' num2str(subj) '/eyetracking/']) '.edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
end
clear resp;


Screen('CloseAll');
ShowCursor;
ListenChar(1);
% catch
%    Screen('CloseAll');
%    ShowCursor;
% end


return