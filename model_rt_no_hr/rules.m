function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

%todo
% hr sensor aansluiten
% als alles werkt: params (b,c) op 0 zetten


% feedback/support =
% in  ~~ pitch ascending
% out ~~ pitch descending


% chest_c in cm = voltage * 10 + 50?
% of logaritmisch?

% ---------------------------------------------------------------
%
%     DOMAIN
%
% ---------------------------------------------------------------

% ! bij waardes uit scenario altijd (t) gebruiken ipv (t+1)
%
% function result = hr( model, trace, parameters, t )
%   hr = 0;
%   result = {t+1, 'hr', hr};
% end
%
% function result = hr_pos( model, trace, parameters, t )
%   % hr_pos = the voltage of the input sensor
%   floor = 1; % param
%
%   global TRAINING
%   if TRAINING
%     global TRAINING_HR
%     pos = TRAINING_HR(t);
%     if pos > floor
%       pos = true;
%     else
%       pos = false;
%     end
%   else
%     pos = false;
%   end
%   result = {t+1, 'hr_pos', pos};
% end

%new
function result = breathing_f( model, trace, parameters, t )
  breathing_f = 0;
  global TRAINING;
  if TRAINING,    breathing_f = 0;  end; %bypass the domain model

  result = {t+1, 'breathing_f', breathing_f};
end

function result = chest_c( model, trace, parameters, t )
  curr_chest_c = 0;
  % new
  global REAL_TIME_INPUT TRAINING TRAINING_BF

  if TRAINING,    curr_chest_c = TRAINING_BF(t)^2 * 10 + 50;  end;
    %todo ^2 weghalen + testen
    %todo deze formula uitleggen of uitwerken
  if REAL_TIME_INPUT
    br = trace(t).breathingvalue.arg{1};
    curr_chest_c = br^4 * 14 + 60;
  end


  % real time graphs
  global CHEST_Y1 PLOT_CHEST1  %CHEST_Y2  RT_CHEST2
  CHEST_Y1(1) = [];
  CHEST_Y1(end+1) = curr_chest_c - 30;
  % CHEST_Y2(end+1) = curr_chest_c;
  refreshdata(PLOT_CHEST1 ); %RT_CHEST1
  % refreshdata(RT_CHEST2);

  result = {t+1, 'chest_c', curr_chest_c};
end
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% ---------------------------------------------------------------
%
%     ANALYSIS
%
% ---------------------------------------------------------------
%
%
%
%

function result = obs_hr( model, trace, parameters, t )
  % for the default domain model
  % hr = trace(t+1).hr.arg{1};
  hr = 75;
  result = {t+1, 'observe', predicate('hr',hr)};
end
% function result = obs_hr_pos( model, trace, parameters, t )
%   % for the training, with real data
%   % hr_pos = the voltage of the input sensor
%   pos = trace(t+1).hr_pos.arg{1};
%   result = {t+1, 'observe', predicate('hr_pos',pos)};
% end

function result = obs_chest_c( model, trace, parameters, t )
  chest_c = trace(t+1).chest_c.arg{1};
  result = {t+1, 'observe', predicate('chest_c',chest_c)};
end

function result = bel_chest_c( model, trace, parameters, t )
  chest_c = l2.getall(trace, t+1, 'observe', predicate('chest_c', NaN)).arg{1}.arg{1};
  % t
  % chest_c
  result = {t+1, 'belief', predicate('chest_c',chest_c)};
end

%new
% function result = bel_hr_pos( model, trace, parameters, t )
%   hr_pos = trace(t+1).hr_pos.arg{1};
%   result = {t+1, 'belief', predicate('hr_pos',hr_pos)};
% end

%new
function result = bel_starting_dir( model, trace, parameters, t )
  if t < 3
    % t
    chest_dir = '1 in';
    % chest_dir
  else
    % t >= 3
    prev2_chest_c = l2.getall(trace,t-1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
    % prev2_chest_c
    prev_chest_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
    % prev_chest_c

    % curr = l2.getall(trace, t+2, 'belief', predicate('chest_c', NaN))
    curr_chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
    % curr_chest_c

    margin  = model.parameters.default.chest_c_margin; % percentage of the range
    max     = model.parameters.default.bel_max_chest_c;
    min     = model.parameters.default.bel_min_chest_c;
    range   = max - min;
    margin  = margin * range; %relative margin

    prev_chest_c = mean([prev2_chest_c prev_chest_c]);

    % if curr_chest_c > prev_chest_c + margin && prev_chest_c > prev2_chest_c + margin
    %   chest_dir = '1 in';
    % elseif curr_chest_c < prev_chest_c - margin && prev_chest_c < prev2_chest_c - margin
    %   chest_dir = '3 out';
    % else
    %   chest_dir = '2 rest';
    % end
    if curr_chest_c > prev_chest_c + margin
      chest_dir = '1 in';
    elseif curr_chest_c < prev_chest_c - margin
    	chest_dir = '3 out';
    else
    	chest_dir = '2 rest';
    end
  end
  result = {t+1, 'belief', predicate('starting_dir',{chest_dir})};
end

function result = bel_breathing_f( model, trace, parameters, t )
  %calculate the believed breathing frequency
  %search backwards; search until enough breathing cycles are found
  n_cycles = model.parameters.default.n_breathing_cycles;
  % max      = model.parameters.default.max_interval_t;
  lo_bf    = model.parameters.default.lowest_bf;
  dt       = model.parameters.default.dt;
  % mode     = model.parameters.default.breathing_dir_mode;

  max_steps = 1/lo_bf / dt;   % 1/0.01 * 1/dt    %new2
  % v will be a vector/list containing 1s and 0s that indicate starting points of breathing cycles.
  v = [];
  a = t;
  mode = 0;
  count = 1;          % matlab starts counting at 1 instead of 0

  % t
  % prev_val3 = l2.getall(trace, a+1, 'belief', predicate('starting_dir', NaN))

  if t == 1
    prev_val = 0;
  else
    prev_val = l2.getall(trace, a+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
  end

  if t <= n_cycles*2  % a minimum of n*2+1 time steps are required (in+out+...+next_in)
    breathing_f = 0;
    % disp('Waiting - collecting data')
    if t==n_cycles*2, disp('  Go  -->');disp('Searching for breathing cycles...');end; %last wait
  else
    while sum(v) <= n_cycles + 1  %the last interval will be cut off
      val = l2.getall(trace, a, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
      % first set the mode (look for inor out)
      if mode == 0
        % Search for the last change in starting_dir (in or out)
        % Only the cycles before this cycle are stored.
        if ~strcmp(val,prev_val)
          if strcmp(val,'1 in')
            mode = 1;
            v = [1];
            t_end = a;              % t_end = t+1 - 1   discard the current t
          elseif strcmp(val,'3 out')
            mode = 3;
            v = [1];
            t_end = a;
          end
        end
      else
        count = count + 1;
        if (strcmp(val, '1 in') && mode==1) || (strcmp(val, '3 out') && mode==3)
          if ~strcmp(val,prev_val)  %store only the start/end of a cycle
            v(count) = 1;
          else
            v(count) = 0;
          end
        else
          v(count) = 0;
        end
      end
      a = a - 1;
      prev_val = val;
      if a <= 1,      break; end; %break at t1
      if count>max_steps,   break; end; %break when max timesteps has been searched
    end % end of while loop

    % Now determine the interval size + number of breathing cycles
    if sum(v) < 2
      breathing_f = 0;
      if t > 6/dt, disp('no cycles found'); end;
    else
      % min n_cycles is not always reached
      % v is a reversed list, thus the last interval is the start (lowest t)
      % the index of the 'start' = the length of the interval
      % e.g. v               = [1 0 0 0 0 1 0 0 0]
      %      corresponding t = [9 8 7 6 5 4 3 2 1]
      %      start or index  = 9 - 4  = 5
      %      used interval   = [1 0 0 0 0]
      %      lenght(interval) = 4
      %      n                = sum(used interval)
      index = length(v) - find(fliplr(v)==1,1);
      v2 = v(1:index);      % n = sum(v2);
      breathing_f = mean(v2) / dt; %new: * 1/dt
    end
  end
  result = {t+1, 'belief', predicate('breathing_f', breathing_f)};
end

%new
function result = bel_hr( model, trace, parameters, t )
  % global TRAINING REAL_TIME_INPUT
  % if ~TRAINING && ~REAL_TIME_INPUT
  %   hr = l2.getall(trace, t+1, 'observe', predicate('hr', NaN)).arg{1}.arg{1};
  % else
  %   % calculate the believed heart rate
  %   % this function functions the same as bel_breathing_f
  %   % search backwards; search until enough cycles are found
  %   n_cycles = model.parameters.default.n_hr_cycles;
  %   dt  = model.parameters.default.dt;
  %   lhr = model.parameters.default.lhr;
  %   % max = 1/dt * model.parameters.default.max_interval_t;
  %   max_steps = 1/lhr / dt;   % 1/0.01 * 1/dt    %new2
  %
  %   % v will be a vector/list containing 1s and 0s that indicate starting points of breathing cycles.
  %   v = [];
  %   a = t;
  %   mode = [];
  %   count = 1;          % matlab starts counting at 1 instead of 0
  %   prev_val = l2.getall(trace, a+1, 'belief', predicate('hr_pos', NaN)).arg{1}.arg{1};
  %
  %   if t <= n_cycles*2  % a minimum of n*2+1 time steps are required (in+out+...+next_in)
  %     hr = 0;
  %     % disp('Waiting - collecting data (hr)')
  %     if t==n_cycles*2, disp('  Go  -->');disp('Searching for heart beats...');end; %last wait
  %   else
  %     while sum(v) <= n_cycles + 1  %the last interval will be cut off
  %       val = l2.getall(trace, a, 'belief', predicate('hr_pos', NaN)).arg{1}.arg{1};
  %       % if val > floor % make the value binary
  %       %   val = 1;
  %       % else
  %       %   val = 0;
  %       % end
  %
  %       % modes:
  %       % if last (newest) observation = 0; mode = 0
  %       % if last (newest) observation > floor; mode = 1
  %       if isempty(mode)
  %         if val ~= prev_val
  %           mode = val; % 1 or 2, true or false
  %           v = [1];
  %           t_end = a;
  %         end
  %       else %mode == 1 or 2
  %         count = count + 1;
  %         if val ~= prev_val
  %           if mode == val
  %             v(count) = 1;
  %           else
  %             v(count) = 0;
  %           end
  %         else
  %           v(count) = 0;
  %         end
  %       end
  %       a = a - 1;
  %       prev_val = val;
  %       if a <= 1,      break; end; %break at t1
  %       if count>max_steps,   break; end; %break when max timesteps has been searched
  %     end % end of while loop
  %
  %     % Now determine the interval size + number of breathing cycles
  %     if sum(v) < 2
  %       hr = 0;
  %       % if t > 6/dt, disp('no beats found'); end;
  %     else
  %       % same as bel bf
  %       % remove partial intervals at the end
  %       index = length(v) - find(fliplr(v)==1,1);
  %       v2 = v(1:index);        % n = sum(v2);
  %       hr = mean(v2) / dt * 60; % in bpm
  %     end
  %   end
  % end
  % global REAL_TIME_INPUT
  % if REAL_TIME_INPUT && hr < lhr, hr = lhr; end;
  hr = 80;
  result = {t+1, 'belief', predicate('hr',hr)};
end

function result = bel_used_chest_range( model, trace, parameters, t )
  %todo: timesteps (review) afhankelijk maken van bf
  %       zodat er naar de laatste 3 breathing cycles gekeken wordt
  bf        = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  min       = model.parameters.default.bel_min_chest_c; % believed min
  max       = model.parameters.default.bel_max_chest_c; % believed max
  dt        = model.parameters.default.dt;
  max_range = max - min;
  highest   = 0;
  lowest    = Inf;

  if bf < 0.01, bf = 0.01; end; %new2, bf heeft al een limiet dus deze kan weg
  time = 1/dt * 1/bf * model.parameters.default.chest_range_cycles;%new2

  if time > t || t < 10 %new2
    relative_range = 1;
  else
    for i=t:-1:t-time
      chest_c = l2.getall(trace, i+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
      if chest_c > highest,   highest = chest_c;          end;
      if chest_c < lowest,    lowest  = chest_c;          end;
    end
    used_range = highest - lowest;
    relative_range = used_range / max_range;
    if relative_range > 1, relative_range = 1;  end;
  end

  % relative_range

  result = {t+1, 'belief', predicate('used_chest_range',relative_range)};
end

function result = bel_br_intensity( model, trace, parameters, t )
  range = l2.getall(trace, t+1, 'belief', predicate('used_chest_range', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('br_intensity',range)};
end

%new
function result = bel_relative_c( model, trace, parameters, t )
  % the previous value of chest_c,
  % independent from the breathing intensity
  %  ?relative to the used chest_range?
  chest_c     = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  min       = model.parameters.default.bel_min_chest_c; % believed min
  max       = model.parameters.default.bel_max_chest_c; % believed max
  range       = max - min;
  A           = range / 2;
  avg_chest_c = min + A;            % 70
  relative_c = (chest_c - avg_chest_c) / A;
    % value is in range [-1,1], except when the range changes between different points in time
  if      relative_c >  1, relative_c =  1;
  elseif  relative_c < -1, relative_c = -1;  end;

  if chest_c > max - 0.5
    relative_c = 1;
  elseif chest_c < min + 0.5
    relative_c = -1;
  end

  result = {t+1, 'belief', predicate('relative_c',relative_c)};
end

% %new
function result = bel_d_bf( model, trace, parameters, t )
  % keep track of whether bf is increasing or decreasing
  % prev_d  = l2.getall(trace, t, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
  prev_bf = l2.getall(trace, t, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  bf      = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  margin = 0.01;
  if bf > prev_bf + margin
    d = 1;
  elseif bf < prev_bf - margin
    d = -1;
  else
    d = 0;
  end
  result = {t+1, 'belief', predicate('d_bf', d)};
end
%new
function result = bel_d_hr( model, trace, parameters, t )
  % qualitative: keep track of whether hr is increasing or decreasing
  % prev_d  = l2.getall(trace, t, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
  prev_hr = l2.getall(trace, t, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  hr      = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  margin = 1;
  if hr > prev_hr + margin
    d = 1;
  elseif hr < prev_hr - margin
    d = -1;
  else
    d = 0;
  end
  result = {t+1, 'belief', predicate('d_hr', d)};
end

% new
function result = bel_anxiety( model, trace, parameters, t )
  prev_anxiety  = l2.getall(trace, t, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  bf            = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  hr            = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  % stable_hr     = l2.getall(trace, t+1, 'belief', predicate('stable_hr', NaN)).arg{1}.arg{1};
  % stable_bf     = l2.getall(trace, t+1, 'belief', predicate('stable_bf', NaN)).arg{1}.arg{1};
  stable_hr = true;
  br_intensity = l2.getall(trace, t+1, 'belief', predicate('br_intensity', NaN)).arg{1}.arg{1};
  % used_chest_range == breathing intensity
  low_int   = model.parameters.default.low_br_int; % b intensity
  dt        = model.parameters.default.dt;
  decay     = model.parameters.default.anxiety_decay;
  floor_bf  = model.parameters.default.floor_bf;
  lhr       = model.parameters.default.lhr;
  a         = model.parameters.default.bf_a;
  b         = model.parameters.default.bf_b;
  c         = model.parameters.default.bf_c;
  pa_time   = 1/dt * model.parameters.default.pa_time;

  % ------------------------------------------------------------------------
  % script to save data (hr + bf) to csv %new2
  % ------------------------------------------------------------------------
  global N TRAINING SAVE_DATA
  if t == N - 1, disp('dt ended with'); disp(dt); end;
  if SAVE_DATA && t == N - 1 %new2
    bb = [];
    hrr = [];
    for i=1:t
        bf3 = l2.getall(trace, i+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
        hr3 = l2.getall(trace, i+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
        bb(i) = bf3;
        hrr(i) = hr3;
    end
    disp('saving traces')
    ra = int2str(rand(1,1)*100);
    name1 = strcat('data/calculated/bb_c2_v',ra,'.csv');
    name2 = strcat('data/calculated/hr_c2_v',ra,'.csv');
    csvwrite(name1,bb);
    csvwrite(name2,hrr);
  end
  % ------------------------------------------------------------------------
  % ------------------------------------------------------------------------

  if t < pa_time
    % the adaption model needs a few timesteps to calculate the parameters for the
    expected_bf = Inf;
  else
    expected_bf = a*hr.^2 + b*hr + c; % f(x) = ax^2 +bx + c; x = x-lhr
  end

  if br_intensity < low_int %&& t>20
    % low intensity
    disp('br_intensity < 0.5 and t=')
    %disp(t)
    anxiety = 10 * (1 - br_intensity);  % 10 * (1 - 0.1) = 90
  % elseif stable_hr && ~stable_bf
  %   % stable hr and unstable bf
  %   disp('stable_hr, unstable bf')
  %   anxiety = 90;
  elseif stable_hr && bf > expected_bf + floor_bf
    % stable hr and higher bf than expected_bf
    disp('stable_hr, bf > expected_bf')
    % expected_bf
    % bf
    anxiety = expected_bf / bf * 100;
  else
    % let the anxiety decay
    anxiety = prev_anxiety * decay;
  end

  if t<50, anxiety = 0; end;
  if anxiety > 100, anxiety = 100; end;
  result = {t+1, 'belief', predicate('anxiety', anxiety)};
end

function result = assessment( model, trace, parameters, t )
  anxiety = l2.getall(trace, t+1, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  floor_a = 0.1;

  if anxiety > floor_a
    assessment = true;
  else
    assessment = false;
  end
  result = {t+1, 'assessment', assessment};
end


%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% ---------------------------------------------------------------
%
%     SUPPORT
%
% ---------------------------------------------------------------
%
%
%
%

function result = des_bf( model, trace, parameters, t )
  % desired breathing_f without influence from anxiety.
  % hr ipv original hr %new2
  hr = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  % h       = model.parameters.default.hr_breathing;
  lhr = model.parameters.default.lhr;
  lo_bf = model.parameters.default.lowest_bf;
  max_bf = model.parameters.default.max_bf;
  a   = model.parameters.default.bf_a;
  b   = model.parameters.default.bf_b;
  c   = 0.8 * model.parameters.default.bf_c; %new2
  % c = c * 0.8;

  % hr = hr_bpm / 60; % in s-1
  bf = a*hr.^2 + b*hr + c;
  if bf < lo_bf
    bf = lo_bf;%new2
  elseif bf > max_bf
    bf = max_bf;%new2
  end

  result = {t+1, 'desire', predicate('breathing_f', bf)};
end

%new
function result = cycle_time( model, trace, parameters, t )
  % we can just give support about breathing in or out (quantitative)
  % not exactly how deep (fast) to breathe (qualitative)
  % on the other hand, we can advice the user to breathe deeper (more intensively)

  % when user should breathe in:  cycle_time = positive
  % when user should breathe out: cycle_time = negative
  % cycle_time is never 0
  prev_ct         = trace(t).cycle_time.arg{1};
  prev_assessment = trace(t).assessment.arg{1};
  assessment      = trace(t+1).assessment.arg{1};
  rel_c           = l2.getall(trace, t+1, 'belief', predicate('relative_c', NaN)).arg{1}.arg{1};
  bf              = l2.getall(trace, t+1, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
  dt              = model.parameters.default.dt;
  max_cycle_time = 1/bf * 1/dt; % 1/0.5 * 1/0.1 = 20

  if ~prev_assessment && assessment
    cycle_time = rel_c * max_cycle_time;
  else

    if prev_ct > 0 % 'in'
      if prev_ct < 0.5 * max_cycle_time % between 0 and 10
        cycle_time = prev_ct + 1;
      else
        cycle_time = -1;
      end
    else % 'out'
      if prev_ct > -0.5 * max_cycle_time % between -10 and 0
        cycle_time = prev_ct - 1;
      else
        cycle_time = 1;
      end
    end

  end
  result = {t+1, 'cycle_time', cycle_time};
end

%new2
function result = des_starting_dir( model, trace, parameters, t )
  cycle_time = trace(t+1).cycle_time.arg{1};
  chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  min       = model.parameters.default.bel_min_chest_c; % believed min
  max       = model.parameters.default.bel_max_chest_c; % believed max
  %todo param adaption gebruiken voor min en max?

  if cycle_time > 0
    if chest_c < max
      dir = '1 in';
    else
      dir = '2 rest';
    end
  else % cycle_time < 0
    if chest_c > min
      dir = '3 out';
    else
      dir = '2 rest';
    end
  end
  result = {t+1, 'desire', predicate('starting_dir', {dir})};
end

function result = support( model, trace, parameters, t )
  starting_dir = l2.getall(trace, t+1, 'desire', predicate('starting_dir', NaN)).arg{1}.arg{1};
  assessment = trace(t+1).assessment.arg{1};
  cycle_time = trace(t+1).cycle_time.arg{1};

  % global PLOT_COLOR PLOT_TXT


  % if assessment
  %   starting_dir = starting_dir;
  %
  %   if strcmp(starting_dir, '1 in')
  %     % PLOT_COLOR(1).FaceColor = [1 0 0];% red
  %     PLOT_TXT = 'data/in.jpg';
  %   elseif strcmp(starting_dir,'3 out')
  %     % PLOT_COLOR(1).FaceColor = [0 0 1];% blue
  %     PLOT_TXT = 'data/uit.jpg';
  %   else
  %     % PLOT_COLOR(1).FaceColor = [1 1 0];% yellow
  %     PLOT_TXT = 'data/rust.jpg';
  %   end
  % else
  %   starting_dir = '4 none';
  %   % PLOT_COLOR(1).FaceColor = [0 0 0];% yellow
  %   PLOT_TXT = 'data/none.jpg';
  % end

  % refreshdata(PLOT_COLOR)
  % imshow(PLOT_TXT);

  global SOUND
  if SOUND && assessment
    dt = model.parameters.default.dt;
    hi_f = 3456; % 8*432    %800 1200
    lo_f = 432; %432 ipv 440
    amp = 0.5;
    % duration = dt - 0.01; %(dt=0.2)
    duration = dt; %(dt=0.2) dt * 2
    Fs = 8192;  % sampling frequency = 8192 2048
      %must be higher than hi_f
    values=0:1/Fs:duration;
    len = length(values);
    decreasing_amp = amp:-1*amp/len:0;
    % a=amp*sin(2*pi * freq * values);
    % sound(a);

    % perfect fifth: 440*1.5= 660 Hz
    % int = 1.5;
    int = 5/4;

    if strcmp(starting_dir,'1 in')
      % f = 800 + (cycle_time*1).^2;
      f = lo_f * int^(cycle_time-1);
      if f > Fs, f = Fs; end; % shannon..
      a2 = decreasing_amp(1:len) .* sin(2*pi* f * values);
      sound(a2,Fs);
    elseif strcmp(starting_dir,'3 out')
      % f = 1200 - (abs(cycle_time)*1).^2;
      % a=amp*sin(2*pi* f * values);
      % sound(a,Fs);
      % f = 1200 - (abs(cycle_time)*1).^2;
      f = hi_f * (1/int)^(abs(cycle_time)-1);
      % a2 = sin(2*pi* f * values);
      % a2 = a2 .* decreasing_amp(1:len);
      if f < 500, f = 500; end; % prevent low bass notes
      a2 = decreasing_amp(1:len) .* sin(2*pi* f * values);
      sound(a2,Fs);
    end
  end

  result = {t+1, 'support', {starting_dir}};
end


%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% --------------------------------------
% --------------------------------------
%
%     PARAMETER ADAPTION
%
% --------------------------------------
% --------------------------------------
%
%

%todo
% gemiddelde waardes gebruiken
%   zoek twee intervallen met een stabiele hr
%   mean(hr) moeten nog steeds ver uit elkaar liggen
%new
% function result = bel_d_hr2( model, trace, parameters, t )
%   % quantitative function, to use with the parameter adaption model
%   prev_hr = l2.getall(trace, t, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
%   hr      = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
%   result = {t+1, 'belief', predicate('d_hr', hr - prev_hr)};
% end
% geen lhr gebruiken, die is overbodig


% andere optie:
% alleen c aanpassen

% function result = adaptions_hr_bf( model, trace, parameters, t )
%   assessment = trace(t+1).assessment.arg{1};
%   skip  = model.parameters.default.pa_skip_n_time_steps;
%   %skip the first 10 timesteps
%   if ~assessment && t>3+skip
%     lhr    = model.parameters.default.lhr;
%     dt     = model.parameters.default.dt;
%     time   = 1/dt * model.parameters.default.pa_time; % the review time
%     margin = 5;
%
%     if time>t,time=t;end;
%     hrs = [];
%     bfs = [];
%     for i=1:3
%       % make random sample of the last x timesteps
%       sample = t+1 - round(time * rand);
%       hr = l2.getall(trace, sample, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
%       bf = l2.getall(trace, sample, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%       hrs(i) = hr;
%       bfs(i) = bf;
%     end
%     if sum(hrs)>lhr*3 && sum(bfs) > 0 && t > 30
%       % hrs
%       % bfs
%       x1 = hrs(1);
%       x2 = hrs(2);
%       y1 = bfs(1);
%       y2 = bfs(2);
%       %todo: hier een interval/average nemen
%
%       if x1 > x2 + margin || x1 < x2 - margin
%         % update the necessary params
%         pa_lim = model.parameters.default.pa_time;
%         a      = model.parameters.default.bf_a;
%         b      = model.parameters.default.bf_b;
%         c      = model.parameters.default.bf_c;
%         s_b    = model.parameters.default.pa_speed_b;
%         s_c    = model.parameters.default.pa_speed_c;
%
%         % a = 1.5505e-04 % when filled in in the 'advanced' formula
%         b2 = (- a*x1^2 + a*x2^2 + y1 - y2)/(x1 - x2);
%         c2 = (a*x1^2*x2 - a*x1*x2^2 + y2*x1 - y1*x2)/(x1 - x2);
%
%         % relative difference of the change
%         rel_d_b = (b2 - b) / b;
%         rel_d_c = (c2 - c) / c;
%
%         % keep diffs in range [-pa_lim,pa_lim]
%         if rel_d_b > pa_lim
%           rel_d_b = pa_lim;
%         elseif rel_d_b < -1*pa_lim
%           rel_d_b = -1*pa_lim;
%         end
%         if rel_d_c > pa_lim
%           rel_d_c = pa_lim;
%         elseif rel_d_c < -1*pa_lim
%           rel_d_c = -1*pa_lim;
%         end
%
%         % apply differences:  p = old p * 1 + (factor * relativ change)
%         b = b * (1 + s_b * rel_d_b);
%         c = c * (1 + s_c * rel_d_c);
%
%         % test the params (bf cannor be < 0)
%         % mini = 30;
%         % maxi = 250;
%         % x = 30:3:250; %linspace
%         % f1 = a*x.^2 + b*x + c;
%         % % f3 = b + 2*a*x;
%         % % index = find(min(f3) == f3)
%         % % lowest = f1(index) % when a>0
%         % if min(f1) > 0
%           %change params
%           model.parameters.default.bf_b = b;
%           model.parameters.default.bf_c = c;
%           % plot/update graph
%           global PLOT_BF HR_AXIS BF_AXIS %BF_A BF_B BF_C
%           BF_AXIS = a*HR_AXIS.^2+b*HR_AXIS+c;
%           % BF_A = BF_A+0.1;
%           % BF_AXIS = BF_A*HR_AXIS;% bf_plot = a*hr_plot.^2+b*hr_plot+c;
%           refreshdata(PLOT_BF)
%         % else
%         %   disp('PA: f1 < 0')
%         % end
%       end
%     end
%   end
%   result = {t+1, 'adaption_1', assessment};
% end

function result = adaptions_chest_c_range( model, trace, parameters, t )
  %adaption min and max chest c
  chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  maxi     = model.parameters.default.bel_max_chest_c;
  mini     = model.parameters.default.bel_min_chest_c;
  change = false;
  if chest_c > maxi
    model.parameters.default.bel_max_chest_c = chest_c;
    change = true;
  end
  if chest_c < mini
    model.parameters.default.bel_min_chest_c = chest_c;
    change = true;
  end
  result = {t+1, 'adaption_2', change};
end

% function result = adaption_dt( model, trace, parameters, t )
%   % change dt when a timestap takes more time
%   dt = model.parameters.default.dt;
%   dt_plus = model.parameters.default.dt_plus;
%   measured_dt = toc;
%   tic;
%
%   global REAL_TIME_INPUT
%   if REAL_TIME_INPUT
%     if dt > measured_dt
%       model.parameters.default.dt = measured_dt + dt_plus;
%       disp('dt++')
%       measured_dt
%     else
%       % measured_dt <= dt
%       % disp(measured_dt - dt)
%       pause(dt - measured_dt);
%     end
%     if mod(t,100) == 0
%       disp('current dt:')
%       measured_dt
%     end
%   end
%   % voor discussie
%   % dt aanpassen leidt tot meetfouten bij o.a. frequenties
%   % dit kan opgelost worden door van dt een concpept te maken
%   % en oude waardes op te slaan in de trace
%   % het bepalen van gemiddelde frequenties kan dan gedaan worden i.c.m. dt
%   result = {t+1, 'adaption_3', false};
% end

function result = adaption_dt( model, trace, parameters, t )
  % change dt when a timestap takes more time
  dt = model.parameters.default.dt;
  dt_plus = model.parameters.default.dt_plus;
  measured_dt = toc;
  tic;

  global REAL_TIME_INPUT LIMIT_DT
  if REAL_TIME_INPUT && t > 20 %new2
    if measured_dt > dt
      if measured_dt > LIMIT_DT + dt_plus, measured_dt = LIMIT_DT; disp('measured_dt > LIMIT_DT!'); end;
      model.parameters.default.dt = measured_dt + dt_plus;
      % disp('dt++')
      % dt
      % measured_dt
    else
      % measured_dt <= dt
      pause(dt - measured_dt);


    end
    if mod(t,100) == 0
      % disp('current dt:')
      % measured_dt
    end
  end
  % voor discussie
  % dt aanpassen leidt tot meetfouten bij o.a. frequenties
  % dit kan opgelost worden door van dt een concpept te maken
  % en oude waardes op te slaan in de trace
  % het bepalen van gemiddelde frequenties kan dan gedaan worden i.c.m. dt
  result = {t+1, 'adaption_3', false};
end


%
%
%
%
%
%
%
%
%
%
%
% ---------------------------------------------------------------
%
%     GRAPHS
%
% ---------------------------------------------------------------
%
%
%
%

% % Breathing_f and heart rate
%
% function result = graph_bel_hr( model, trace, parameters, t )
%   bf = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_bel_hr', bf};
% end
function result = graph_bel_breathing_f( model, trace, parameters, t )
  bf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_breathing_f', bf};
end
% function result = graph_breathing_f_error( model, trace, parameters, t )
%   bf  = l2.getall(trace, t+1, 'breathing_f', NaN).arg{1};
%   bbf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_breathing_f_error', bbf-bf};
% end
% % function result = graph_d_hr( model, trace, parameters, t )
% %   x = l2.getall(trace, t+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
% %   result = {t+1, 'graph_d_hr', x};
% % end
% % function result = graph_d_bf( model, trace, parameters, t )
% %   x = l2.getall(trace, t+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
% %   result = {t+1, 'graph_d_bf', x};
% % end
%
% % Starting direction, physical state, original_hr, anxiety
%
function result = graph_bel_starting_dir( model, trace, parameters, t )
  if t == 1
    chest_pos = '1 in';
  else
    chest_pos = l2.getall(trace, t+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
  end
  result = {t+1, 'graph_bel_starting_dir', {chest_pos}};
end
function result = graph_bel_anxiety( model, trace, parameters, t )
  h = l2.getall(trace, t+1, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_anxiety', h};
end
% % function result = graph_bel_prev_ps( model, trace, parameters, t )
% %   ps = l2.getall(trace, t+1, 'belief', predicate('ps', NaN)).arg{1}.arg{1};
% %   result = {t+1, 'graph_bel_prev_ps', ps};
% % end
% % function result = graph_original_hr( model, trace, parameters, t )
% %   h = l2.getall(trace, t+1, 'belief', predicate('original_hr', NaN)).arg{1}.arg{1};
% %   result = {t+1, 'graph_original_hr', h};
% % end
% function result = graph_bel_used_chest_range( model, trace, parameters, t )
%   x = l2.getall(trace, t+1, 'belief', predicate('used_chest_range', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_bel_used_chest_range', x};
% end
%
% function result = graph_stable_hr( model, trace, parameters, t )
%   x = l2.getall(trace, t+1, 'belief', predicate('stable_hr', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_stable_hr', x};
% end
% function result = graph_stable_bf( model, trace, parameters, t )
%   x = l2.getall(trace, t+1, 'belief', predicate('stable_bf', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_stable_bf', x};
% end
%
% % function result = graph_bel_phi( model, trace, parameters, t )
% %   x   = l2.getall(trace, t+1, 'belief', predicate('phase_shift', NaN)).arg{1}.arg{1};
% %   result = {t+1, 'graph_bel_phi', x};
% % end
%
% % Support model

function result = graph_des_bf( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_des_breathing_f', x};
end
function result = graph_breathing_f_diff( model, trace, parameters, t )
  bel_bf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  des_bf = l2.getall(trace, t+1, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_breathing_f_diff', bel_bf - des_bf};
end
function result = graph_des_starting_dir( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'desire', predicate('starting_dir', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_des_starting_dir', x};
end
