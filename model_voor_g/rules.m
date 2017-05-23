function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end


% chest_c in cm = voltage * 10 + 50?
% of logaritmisch?


%todo:
% dt 'vastzetten' op e.g. dt = 0.19

% of regel bij min voltage
% if voltage < min_voltage
%   min_voltage = voltage - 0.01
% idem voor max
% op deze manier kan je instellen op 0 en gaan ze vanzelf goed staan
% mits de user 1x diep in en uitademt



% start t=0 with tic
% before going to t=1, check tic. Two options
% (1) if tac < dt, pause( dt - tac ); else, disp('taking too much time') (or raise dt)
% (2) dt(t+1) = tac

% ---------------------------------------------------------------
%
%     DOMAIN
%
% ---------------------------------------------------------------

% ! bij waardes uit scenario altijd (t) gebruiken ipv (t+1)

%todo: geen prev_... concepten gebruiken

function result = mind( model, trace, parameters, t )
  % dir = trace(t).support.arg{1};
  dir = '4 none';
  result = {t+1, 'mind', {dir}};
end

function result = anxiety_regulation( model, trace, parameters, t )
  % 0 = bad, 1 = perfect
  breathing_f = trace(t).breathing_f.arg{1};
  disfac      = model.parameters.default.disfac_a_reg;
  b           = model.parameters.default.breathing_anxiety;
  % lo_bf = alwas zero
  hi_bf       = model.parameters.default.high_bf;
  relative_bf = breathing_f / hi_bf;
  anxiety_reg = disfac - b * relative_bf;
  if anxiety_reg <= 0, anxiety_reg = 0; disp('lowest anxiety regulation reached');  end;
  result = {t+1, 'anxiety_regulation', anxiety_reg};
end

function result = anxiety( model, trace, parameters, t )
  % percentage of full capacity
  sitfac        = trace(t).sitfac.arg{1}; %gives errors (t+1 of t)
  prev_anxiety  = trace(t).anxiety.arg{1};
  regulation    = trace(t+1).anxiety_regulation.arg{1};
  decay         = model.parameters.default.anxiety_decay;
  disfac        = model.parameters.default.disfac;
  % Inf waardes geven errors... wachten op bugfix van linford

  anxiety = (disfac * sitfac + prev_anxiety * decay) * (1 - regulation);
  result = {t+1, 'anxiety', anxiety};
end

function result = hr_var( model, trace, parameters, t )
  % hr_var = the amount of bpm that will be 'added' or 'subtracted'
  anxiety  = trace(t+1).anxiety.arg{1};
  d        = model.parameters.default.disfac_hr_var;
  a        = model.parameters.default.anxiety_hr_var;
  % influence from ps..... = 0

  hr_var = d * rand + a * anxiety * rand;
  result = {t+1, 'hr_var', hr_var};
end

function result = hr( model, trace, parameters, t )
  % in bpm, between 0 and 250
  % output = hr * dt
  ps      = trace(t).ps.arg{1}; % t+1 doesn't work with this syntax + scenario values
  var     = trace(t+1).hr_var.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  a       = model.parameters.default.anxiety_hr;
  bhr     = model.parameters.default.bhr;
  lhr     = model.parameters.default.lhr;

  hr = (bhr + ps) + (a * anxiety) + var;
  if hr < lhr, hr = lhr; disp('lhr reached!'); end;

  global TRAINING;
  if TRAINING,    hr = 0;  end; %bypass the domain model
  result = {t+1, 'hr', hr};
end

function result = hr_pos( model, trace, parameters, t )
  % hr_pos = the voltage of the input sensor
  floor = 1; % param

  global TRAINING;
  if TRAINING
    global TRAINING_HR
    pos = TRAINING_HR(t);
    if pos > floor
      pos = true;
    else
      pos = false;
    end
  else
    pos = false;
  end
  result = {t+1, 'hr_pos', pos};
end

%new
function result = breathing_f( model, trace, parameters, t )
  % between 0 and 4 (params.high_bf)
  % hr does not have to be translated to bpm
  hr      = trace(t+1).hr.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  a2      = model.parameters.default.anxiety_bf;
  %h       = model.parameters.default.hr_breathing;
  % h2      = model.parameters.default.hr_breathing_exp;
  lhr     = model.parameters.default.lhr;
  a       = model.parameters.default.default_a;
  b       = model.parameters.default.default_b;
  c       = model.parameters.default.default_c;
  % hr = hr_bpm / 60; % in s-1
  % breathing_f = a*hr.^2 + b*hr + c + (a2 * anxiety);
  breathing_f = a*(hr-lhr)^2 + b*(hr-lhr) + c + (a2 * anxiety);

  global TRAINING;
  if TRAINING,    breathing_f = 0;  end; %bypass the domain model

  result = {t+1, 'breathing_f', breathing_f};
end

%new
function result = br_intensity( model, trace, parameters, t )
  % range = [0;1]
  dispos = trace(t).dispos_br_i.arg{1};  %t+1 doesn't work for scenario predicates
  anxiety = trace(t+1).anxiety.arg{1};
  a = model.parameters.default.anxiety_br_i;

  br_intensity = dispos - a * anxiety;      % e.g. 0.9 - 0.5 = 0.4 = 40%
  if br_intensity < 0, br_intensity = 0;  end;
  result = {t+1, 'br_intensity', br_intensity};
end

%new
function result = used_chest_range( model, trace, parameters, t )
  % range = [0;1]
  % if higher than 1; subject is breathing more intensely than healthy
  br_intensity = trace(t+1).br_intensity.arg{1};  %t+1 doesn't work for scenario predicates

  max = model.parameters.default.max_chest_range;
  used_chest_range = max * br_intensity;
  result = {t+1, 'used_chest_range', used_chest_range};
end

function result = chest_c( model, trace, parameters, t )

  global TRAINING;
  global TRAINING_BF;
  if TRAINING,    curr_chest_c = TRAINING_BF(t)^2 * 10 + 50;  end;
    %todo deze formula uitleggen of uitwerken



  result = {t+1, 'chest_c', curr_chest_c};
end





%new relative_chest_c ipv prev...
function result = relative_c( model, trace, parameters, t )
  %the previous value of chest_c, relative to the preferred/used chest_range
  chest_c     = trace(t+1).chest_c.arg{1};            % 60-80
  range       = trace(t+1).used_chest_range.arg{1}; % 19
  min         = model.parameters.default.min_chest_c;
  max         = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2);            % 70
  A = range / 2;

  relative_c = (chest_c - avg_chest_c) / A;
    % value is in range [-1,1], except when the range changes between different points in time
  if      relative_c >  1, relative_c =  1;
  elseif  relative_c < -1, relative_c = -1;  end;
  result = {t+1, 'relative_c', relative_c};
end
%new

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
  hr = trace(t+1).hr.arg{1};
  result = {t+1, 'observe', predicate('hr',hr)};
end
function result = obs_hr_pos( model, trace, parameters, t )
  % for the training, with real data
  % hr_pos = the voltage of the input sensor
  pos = trace(t+1).hr_pos.arg{1};
  result = {t+1, 'observe', predicate('hr_pos',pos)};
end

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
function result = bel_hr_pos( model, trace, parameters, t )
  hr_pos = trace(t+1).hr_pos.arg{1};
  result = {t+1, 'belief', predicate('hr_pos',hr_pos)};
end

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
  max      = model.parameters.default.max_interval_t;
  dt       = model.parameters.default.dt;
  % mode     = model.parameters.default.breathing_dir_mode;

  max = max / dt;   % t=1,60x     t=0.1,600x
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
    disp('Waiting - collecting data')
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
      if count>max,   break; end; %break when max timesteps has been searched
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
  global TRAINING
  if ~TRAINING
    hr = l2.getall(trace, t+1, 'observe', predicate('hr', NaN)).arg{1}.arg{1};
  else
    % calculate the believed heart rate
    % this function functions the same as bel_breathing_f
    % search backwards; search until enough cycles are found
    n_cycles = model.parameters.default.n_hr_cycles;
    dt  = model.parameters.default.dt;
    max = 1/dt * model.parameters.default.max_interval_t; % same as bel bf

    % v will be a vector/list containing 1s and 0s that indicate starting points of breathing cycles.
    v = [];
    a = t;
    mode = [];
    count = 1;          % matlab starts counting at 1 instead of 0
    prev_val = l2.getall(trace, a+1, 'belief', predicate('hr_pos', NaN)).arg{1}.arg{1};

    if t <= n_cycles*2  % a minimum of n*2+1 time steps are required (in+out+...+next_in)
      hr = 0;
      disp('Waiting - collecting data (hr)')
      if t==n_cycles*2, disp('  Go  -->');disp('Searching for heart beats...');end; %last wait
    else
      while sum(v) <= n_cycles + 1  %the last interval will be cut off
        val = l2.getall(trace, a, 'belief', predicate('hr_pos', NaN)).arg{1}.arg{1};
        % if val > floor % make the value binary
        %   val = 1;
        % else
        %   val = 0;
        % end

        % modes:
        % if last (newest) observation = 0; mode = 0
        % if last (newest) observation > floor; mode = 1
        if isempty(mode)
          if val ~= prev_val
            mode = val; % 1 or 2, true or false
            v = [1];
            t_end = a;
          end
        else %mode == 1 or 2
          count = count + 1;
          if val ~= prev_val
            if mode == val
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
        if count>max,   break; end; %break when max timesteps has been searched
      end % end of while loop

      % Now determine the interval size + number of breathing cycles
      if sum(v) < 2
        hr = 0;
        if t > 6/dt, disp('no beats found'); end;
      else
        % same as bel bf
        % remove partial intervals at the end
        index = length(v) - find(fliplr(v)==1,1);
        v2 = v(1:index);        % n = sum(v2);
        hr = mean(v2) / dt * 60; % in bpm
      end
    end

  end
  result = {t+1, 'belief', predicate('hr',hr)};
end







% new
function result = bel_anxiety( model, trace, parameters, t )
  prev_anxiety  = l2.getall(trace, t, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  bf            = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  hr            = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  % stable_hr     = l2.getall(trace, t+1, 'belief', predicate('stable_hr', NaN)).arg{1}.arg{1};
  % stable_bf     = l2.getall(trace, t+1, 'belief', predicate('stable_bf', NaN)).arg{1}.arg{1};
  % br_intensity = l2.getall(trace, t+1, 'belief', predicate('br_intensity', NaN)).arg{1}.arg{1};
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

  % script om data op te slaan ---------------------------------------------
  global N TRAINING
  if t == N - 1 && TRAINING
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
    %   xlswrite('data/hee.xls',bb);
  end
  % ------------------------------------------------------------------------


  % if t < pa_time
  %   % the adaption model needs a few timesteps to calculate the parameters for the
  %   expected_bf = Inf;
  % else
  %   expected_bf = a*(hr-lhr).^2 + b*(hr-lhr) + c; % f(x) = ax^2 +bx + c; x = x-lhr
  % end

  % if br_intensity < low_int
  %   % low intensity
  %   disp('br_intensity < 0.5')
  %   anxiety = 10 * (1 - br_intensity);  % 10 * (1 - 0.1) = 90
  % elseif stable_hr && ~stable_bf
  %   % stable hr and unstable bf
  %   disp('stable_hr, unstable bf')
  %   anxiety = 90;
  % elseif stable_hr && bf > expected_bf + floor_bf
  %   % stable hr and higher bf than expected_bf
  %   disp('stable_hr, bf > expected_bf')
  %   % expected_bf
  %   % bf
  %   anxiety = expected_bf / bf * 100;
  % else
  %   % let the anxiety decay
  %   anxiety = prev_anxiety * decay;
  % end

  % if anxiety > 100, anxiety = 100; end;
  anxiety = 10;
  result = {t+1, 'belief', predicate('anxiety', anxiety)};
end

function result = assessment( model, trace, parameters, t )
  anxiety = l2.getall(trace, t+1, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  floor_a = 0.0001;
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
% ---------------------------------------------------------------
%
%     PARAMETER ADAPTION
%
% ---------------------------------------------------------------
%
%
%
%


function result = adaptions_chest_c_range( model, trace, parameters, t )
  %adaption min and max chest c
  chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  max     = model.parameters.default.bel_max_chest_c;
  min     = model.parameters.default.bel_min_chest_c;
  change = false;
  if chest_c > max
    model.parameters.default.bel_max_chest_c = chest_c;
    change = true;
  end
  if chest_c < min
    model.parameters.default.bel_min_chest_c = chest_c;
    change = true;
  end
  % max
  % min
  result = {t+1, 'adaption_2', change};
end


% if TRAINING
%   dt = fixed
% else
%   if dt2 > dt
%     dt++
% wait until dt
% go


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

% Breathing_f and heart rate

function result = graph_bel_hr( model, trace, parameters, t )
  bf = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_hr', bf};
end
function result = graph_bel_breathing_f( model, trace, parameters, t )
  bf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_breathing_f', bf};
end
function result = graph_breathing_f_error( model, trace, parameters, t )
  bf  = l2.getall(trace, t+1, 'breathing_f', NaN).arg{1};
  bbf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_breathing_f_error', bbf-bf};
end
% function result = graph_d_hr( model, trace, parameters, t )
%   x = l2.getall(trace, t+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_d_hr', x};
% end
% function result = graph_d_bf( model, trace, parameters, t )
%   x = l2.getall(trace, t+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_d_bf', x};
% end

% Starting direction, physical state, original_hr, anxiety
