function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end


% chest_c in cm = voltage * 10 + 50?
% of logaritmisch?

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
  % hr does not have to be translated to bpm because the formula with ax bx c is used
  hr      = trace(t+1).hr.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  a2      = model.parameters.default.anxiety_bf;
  lhr     = model.parameters.default.lhr;
  a       = model.parameters.default.default_a;
  b       = model.parameters.default.default_b;
  c       = model.parameters.default.default_c;
  % hr = hr_bpm / 60; % in s-1
  % breathing_f = a*hr.^2 + b*hr + c + (a2 * anxiety);
  breathing_f = a*(hr)^2 + b*(hr) + c + (a2 * anxiety);

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
  chest_c = trace(t).chest_c.arg{1};
  f       = trace(t).breathing_f.arg{1};
  range   = trace(t).used_chest_range.arg{1};
  phi     = trace(t).phase_shift.arg{1};
  dt          = model.parameters.default.dt;
  min         = model.parameters.default.min_chest_c;
  max         = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2); % 80
  % f = dt * f; % breathing_f is already inheritably adapted to dt,
                % but if dt is changed during simulation, it could produce weird behaviour.
  A = range / 2;
  % y(t) = A sin(2 pi f t dt + phi)
  curr_chest_c = A * sin(2*pi* f * dt + phi) + avg_chest_c;

  %new
  % assume that user always listens to support
  dir = trace(t+1).mind.arg{1}; %new
  if strcmp(dir,'1 in')
    if chest_c < max
      curr_chest_c = chest_c + (max - chest_c) / 2;
      %note that max will never will reached
    else
      curr_chest_c = max;
    end
  elseif strcmp(dir,'3 out')
    if chest_c > min
      curr_chest_c = chest_c - (chest_c - min) / 2;
    else
      curr_chest_c = min;
    end
  elseif strcmp(dir,'2 rest')
    curr_chest_c = chest_c;
  end

  % new
  global TRAINING;
  global TRAINING_BF;
  if TRAINING,    curr_chest_c = TRAINING_BF(t)^2 * 10 + 50;  end;

  % real time graphs
  global CHEST_Y1 PLOT_CHEST1  %CHEST_Y2  RT_CHEST2
  CHEST_Y1(1) = [];
  CHEST_Y1(end+1) = curr_chest_c - 50;
  % CHEST_Y2(end+1) = curr_chest_c;
  refreshdata(PLOT_CHEST1 ); %RT_CHEST1
  % refreshdata(RT_CHEST2);

  result = {t+1, 'chest_c', curr_chest_c};
end

function result = starting_dir( model, trace, parameters, t )
  % in the domain model the starting_dir is calculated instead of observed
  %
  % calculate the direction of the breathing cycle on a point t
  % this is done by comparing x(t) and x(t+dt)
  % to be precise, dt has to be as small as possible (lim->0) thus a new value for dt will be used.
  % (t-dt) and (t+dt) is more realistic, but this can prevent the cycle from completing in certain instance.
  t1 = 1;
  t2 = 1.0001;
  phi   = trace(t).phase_shift.arg{1};
  f     = trace(t).breathing_f.arg{1};
  dt    = model.parameters.default.dt;
  % x1 = A * sin(2*pi* f * t1 + phi) + avg_chest_c;
  % x2 = A * sin(2*pi* f * t2 + phi) + avg_chest_c;
  % A, avg_chest_c can be disregarded
  x1 = sin(2*pi*f*dt* t1 + phi);
  x2 = sin(2*pi*f*dt* t2 + phi);

  if x1 < x2
    curr_dir2 = '1 in';
  else
    curr_dir2 = '3 out';
  end

  result = {t+1, 'starting_dir', {curr_dir2}};
end

function result = performance( model, trace, parameters, t )
  % percentage of full capacity
  anxiety  = trace(t+1).anxiety.arg{1};
  performance = 100 - anxiety;
  if performance < 0, performance = 0;  end;
    % disp('performance')
  result = {t+1, 'performance', {performance}};
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

%new prev...
function result = phase_shift( model, trace, parameters, t )
  %phase_shift of chest_c-cycle on t-1,
  % value will be used in (t-1) to calculate chest_c
  % to calculate the phase, we need the new used_chest_range

  prev_starting_dir = trace(t+1).starting_dir.arg{1};
  prev_relative_c   = trace(t+1).relative_c.arg{1};
  range             = trace(t+1).used_chest_range.arg{1};

  % Formula for a basic sin wave: x(t) = A * sin(2*pi* f * t) + d;
  % However, f can change during a cycle
  % Adapt the formula so that f can be changed:
  %   Let t = 1 and use +phi to shift the wave (with +phi in range[0,2pi])
  %   phi has to be calculate manually because the range is also variable
  %   (phi = asin( x(t)/prev_A) -2*pi*f*t) can not be used

  % calculate the phase shift (phi)
  %  2pi = full shift
  %  1pi = half shift (inverted direction)
  % .5pi = heighest point

  % val = relative_c * 0.5 * pi;
  % omschrijven: use asin to convert the [-1,1] to a 'sin-shape'
  % want sin(t*0.5pi) is niet evenredig met t[0:1]
  % val = asin(relative_c)/pi*2    * 0.5 * pi
  %     Dit  voorkomt ook het vastlopen van de grafiek op relative_c = 1 of -1
  val = asin(prev_relative_c);

  if strcmp(prev_starting_dir, '1 in')
    % phi = relative_c * 0.5 * pi;
    phi = val;
  else
    %   % prev_starting_dir = '3 out'
    % phi = relative_c * -0.5*pi + pi;
    phi = val  * -1 + pi;
    % eigenlijk:
    % if relative_c > 0       % make pos relative_c negative
    %   phi = val * -1 + pi;
    % else                    % make negative relative_c positive
    %   phi = val * -1 + pi;
    % end
  end
  % relative_c
  % phi
  result = {t+1, 'phase_shift', phi}; % + 0.25*pi
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
      % disp('Waiting - collecting data (hr)')
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

function result = bel_used_chest_range( model, trace, parameters, t )
  %todo: timesteps (review) afhankelijk maken van bf
  %       zodat er naar de laatste 3 breathing cycles gekeken wordt
  min       = model.parameters.default.bel_min_chest_c; % believed min
  max       = model.parameters.default.bel_max_chest_c; % believed max
  dt        = model.parameters.default.dt;
  time      = 1/dt * model.parameters.default.chest_range_time;
  max_range = max - min;
  highest   = 0;
  lowest    = Inf;

  if time > t
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

%new
%todo: deze zijn nog altijd 0
function result = bel_stable_bf( model, trace, parameters, t )
  % check if bf has both increased and decreased in a certain interval
  bf        = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  dt        = model.parameters.default.dt;
  n_cycles  = 0.2; % number of cycles to search
  time_per_cycle = bf / dt;
  time = n_cycles * time_per_cycle;
  if time>t, time=t; end;
  positives = [0];
  negatives = [0];
  for i=t:-1:t - time
    val = l2.getall(trace, i+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
    if val == 1
      positives(end) = 1;
    elseif val == -1
      negatives(end) = 1;
    end
  end

  if length(positives) > 1 && length(negatives) > 1
    stable = false;
  else
    stable = true;
  end
  result = {t+1, 'belief', predicate('stable_bf', stable)};
end
%new
function result = bel_stable_hr( model, trace, parameters, t )
  % check if hr has both increased and decreased in a certain interval
  % hr        = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  dt   = model.parameters.default.dt;
  time = 1/dt * model.parameters.default.dt; % number of cycles to search
  if time>t, time=t; end;
  positives = [0];
  negatives = [0];
  for i=t:-1:t - time
    val = l2.getall(trace, i+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
    if val == 1
      positives(end) = 1;
    elseif val == -1
      negatives(end) = 1;
    end
  end

  if length(positives) > 1 && length(negatives) > 1
    stable = false;
  else
    stable = true;
  end
  result = {t+1, 'belief', predicate('stable_hr', stable)};
end


% new
function result = bel_anxiety( model, trace, parameters, t )
  prev_anxiety  = l2.getall(trace, t, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  bf            = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  hr            = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  stable_hr     = l2.getall(trace, t+1, 'belief', predicate('stable_hr', NaN)).arg{1}.arg{1};
  stable_bf     = l2.getall(trace, t+1, 'belief', predicate('stable_bf', NaN)).arg{1}.arg{1};
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

  % script to save data as .csv file ---------------------------------------------
  % global N TRAINING
  % if t == N - 1 && TRAINING
  %   bb = [];
  %   hrr = [];
  %   for i=1:t
  %       bf3 = l2.getall(trace, i+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  %       hr3 = l2.getall(trace, i+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  %       bb(i) = bf3;
  %       hrr(i) = hr3;
  %   end
  %   disp('saving traces')
  %   ra = int2str(rand(1,1)*100);
  %   name1 = strcat('data/calculated/bb_c1_v',ra,'.csv');
  %   name2 = strcat('data/calculated/hr_c1_v',ra,'.csv');
  %   csvwrite(name1,bb);
  %   csvwrite(name2,hrr);
  %   %   xlswrite('data/hee.xls',bb);
  % end
  % ------------------------------------------------------------------------


  if t < pa_time
    % the adaption model needs a few timesteps to calculate the parameters for the
    expected_bf = Inf;
  else
    expected_bf = a*hr.^2 + b*hr + c; % f(x) = ax^2 +bx + c
  end

  if br_intensity < low_int
    % low intensity
    disp('br_intensity < 0.5')
    anxiety = 10 * (1 - br_intensity);  % 10 * (1 - 0.1) = 90
  elseif stable_hr && ~stable_bf
    % stable hr and unstable bf
    disp('stable_hr, unstable bf')
    anxiety = 90;
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

  if anxiety > 100, anxiety = 100; end;
  result = {t+1, 'belief', predicate('anxiety', anxiety)};
end

function result = bel_ps( model, trace, parameters, t )
  anxiety = l2.getall(trace, t+1, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  hr      = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  a       = model.parameters.default.anxiety_hr;
  bhr     = model.parameters.default.bhr;
  % omschrijven: voor in appendix
  % hr = (bhr * ps) + (a * anxiety);
  % hr = (bhr * ps) + (a * anxiety)
  % hr - a * anxiety = bhr * ps
  % ps = (hr - a * anxiety) / bhr
  ps = (hr - a * anxiety) / bhr;  % optional: + param: deviation of ps?
  result = {t+1, 'belief', predicate('ps', ps)};
end

function result = bel_original_hr( model, trace, parameters, t )
  %the value of the heart rate without the influence from anxiety
  ps  = l2.getall(trace, t+1, 'belief', predicate('ps', NaN)).arg{1}.arg{1};
  bhr = model.parameters.default.bhr; %todo kan dit(2)?
  % hr = (bhr * ps) + (a * anxiety)
  hr = bhr * ps;
  result = {t+1, 'belief', predicate('original_hr', hr)};
end

function result = assessment( model, trace, parameters, t )
  anxiety = l2.getall(trace, t+1, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  floor_a = 0.0001;

  %demo
  % if t < 100
  % assessment = true

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
  % hr ipv original hr
  hr = l2.getall(trace, t+1, 'belief', predicate('original_hr', NaN)).arg{1}.arg{1};
  % h       = model.parameters.default.hr_breathing;
  lhr = model.parameters.default.lhr;
  a   = model.parameters.default.bf_a;
  b   = model.parameters.default.bf_b;
  c   = model.parameters.default.bf_c;

  % hr = hr_bpm / 60; % in s-1
  breathing_f = a*hr.^2 + b*hr + c;

  %testing:   this should determine the final bf of the user
  % breathing_f = 0.4;
  result = {t+1, 'desire', predicate('breathing_f', breathing_f)};
end



% oude des_starting_dir werkt niet:
%   support is afhankelijk van chest c, maar stuurt de breathing intensity
%   (snelheid van het ademen zelf) niet aan.
%
%

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

  global PLOT_COLOR

  % assessment = true; %testing

  if assessment
    starting_dir = starting_dir;

    if strcmp(starting_dir, '1 in')
      PLOT_COLOR(1).FaceColor = [1 0 0];% red
    elseif strcmp(starting_dir,'3 out')
      PLOT_COLOR(1).FaceColor = [0 0 1];% blue
    else
      PLOT_COLOR(1).FaceColor = [1 1 0];% yellow
      % h(1).FaceColor = [0 0.5 0];% green
    end
  else
    starting_dir = '4 none';
    PLOT_COLOR(1).FaceColor = [0 0 0];% yellow
  end

  refreshdata(PLOT_COLOR)

  %------------------------------------------------------------------
  %-------- Real time synthesis - can be igored in the rules --------
  %------------------------------------------------------------------
  global SOUND
  if SOUND
    dt = model.parameters.default.dt;
    hi_f = 3456; % 8*432    %800 1200
    lo_f = 432; %432 ipv 440
    amp = 0.5;
    duration = dt - 0.05; %(dt=0.2)
    Fs = 8192;  % sampling frequency = 8192 2048
      %must be higher than hi_f
    values = 0:1/Fs:duration;
    len = length(values);
    decreasing_amp = amp:-1*amp/len:0;
    % perfect fifth: 440*1.5= 660 Hz
    % int = 1.5;
    int = 5/4;

    if strcmp(starting_dir,'1 in')
      % play notes with an increasing pitch

      % f = 800 + (cycle_time*1).^2;
      f = lo_f * int^(cycle_time-1);
      if f > Fs, f = Fs; end; % shannon..
      a2 = decreasing_amp(1:len) .* sin(2*pi* f * values);
      sound(a2,Fs);

    elseif strcmp(starting_dir,'3 out')
      % play notes with a decreasing pitch

      % a=amp*sin(2*pi* f * values);
      f = hi_f * (1/int)^(abs(cycle_time)-1);
      if f < 500, f = 500; end; % prevent low bass notes
      a2 = decreasing_amp(1:len) .* sin(2*pi* f * values);
      sound(a2,Fs);

    end
  end
  %------------------------------------------------------------------
  %------------------------------------------------------------------
  %------------------------------------------------------------------

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
% ---------------------------------------------------------------
%
%     PARAMETER ADAPTION
%
% ---------------------------------------------------------------
%
%
%
%

%todo
% gemiddelde waardes gebruiken
%   zoek twee intervallen met een stabiele hr
%   mean(hr) moeten nog steeds ver uit elkaar liggen

function result = adaptions_hr_bf( model, trace, parameters, t )
  assessment = trace(t+1).assessment.arg{1};
  skip  = model.parameters.default.pa_skip_n_time_steps;
  %skip the first 10 timesteps and skip a minimum of 3 timesteps
  if ~assessment && t>3+skip
    lhr    = model.parameters.default.lhr;
    dt     = model.parameters.default.dt;
    time   = 1/dt * model.parameters.default.pa_time; % the review time
    margin = 5;

    if time>t,time=t;end;
    hrs = [];
    bfs = [];
    for i=1:3
      % make random sample of the last x timesteps
      sample = t+1 - round(time * rand);
      hr = l2.getall(trace, sample, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
      bf = l2.getall(trace, sample, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
      hrs(i) = hr;
      bfs(i) = bf;
    end
    if sum(hrs)>lhr*3 && sum(bfs) > 0 && t > 30

      % disp(length(hrs))
      x1 = hrs(1);    y1 = bfs(1);
      x2 = hrs(2);    y2 = bfs(2);
      %x3 and %y3 are unused

      if x1 > x2 + margin || x1 < x2 - margin
        % update the necessary params
        pa_lim = model.parameters.default.pa_time;
        a      = model.parameters.default.bf_a;
        b      = model.parameters.default.bf_b;
        c      = model.parameters.default.bf_c;
        s_b    = model.parameters.default.pa_speed_b;
        s_c    = model.parameters.default.pa_speed_c;

        % a = 1.5505e-04 % when filled in in the 'advanced' formula
        b2 = (- a*x1^2 + a*x2^2 + y1 - y2)/(x1 - x2);
        c2 = (a*x1^2*x2 - a*x1*x2^2 + y2*x1 - y1*x2)/(x1 - x2);

        % relative difference of the change
        rel_d_b = (b2 - b) / b;
        rel_d_c = (c2 - c) / c;

        % keep diffs in range [-pa_lim,pa_lim]
        if rel_d_b > pa_lim
          rel_d_b = pa_lim;
        elseif rel_d_b < -1*pa_lim
          rel_d_b = -1*pa_lim;
        end
        if rel_d_c > pa_lim
          rel_d_c = pa_lim;
        elseif rel_d_c < -1*pa_lim
          rel_d_c = -1*pa_lim;
        end

        % apply differences:  p = old p * 1 + (factor * relative change)
        b = b * (1 + s_b * rel_d_b);
        c = c * (1 + s_c * rel_d_c);

        %change params
        model.parameters.default.bf_b = b;
        model.parameters.default.bf_c = c;

        % plot/update graph
        global PLOT_BF HR_AXIS BF_AXIS %BF_A BF_B BF_C
        BF_AXIS = a*HR_AXIS.^2+b*HR_AXIS+c;
        % BF_A = BF_A+0.1;
        % BF_AXIS = BF_A*HR_AXIS;% bf_plot = a*hr_plot.^2+b*hr_plot+c;
        refreshdata(PLOT_BF)
      end
    end
  end
  result = {t+1, 'adaption_1', assessment};
end

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

function result = adaption_dt( model, trace, parameters, t )
  % change dt when a timestap takes more time
  dt = model.parameters.default.dt;
  dt_plus = model.parameters.default.dt_plus;
  measured_dt = toc;
  tic;

  global REAL_TIME_INPUT
  if REAL_TIME_INPUT
    if dt > measured_dt
      model.parameters.default.dt = measured_dt + dt_plus;
      disp('dt++')
      measured_dt
    else
      % measured_dt <= dt
      % disp(measured_dt - dt)
      pause(dt - measured_dt);
    end
    if mod(t,100) == 0
      disp('current dt:')
      measured_dt
    end
  end
  % voor discussie
  % dt aanpassen leidt tot meetfouten bij o.a. frequenties
  % dit kan opgelost worden door van dt een concpept te maken
  % en oude waardes op te slaan in de trace
  % het bepalen van gemiddelde frequenties kan dan gedaan worden i.c.m. dt
  result = {t+1, 'adaption_3', false};
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
function result = graph_bel_prev_ps( model, trace, parameters, t )
  ps = l2.getall(trace, t+1, 'belief', predicate('ps', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_prev_ps', ps};
end
function result = graph_original_hr( model, trace, parameters, t )
  h = l2.getall(trace, t+1, 'belief', predicate('original_hr', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_original_hr', h};
end
function result = graph_bel_used_chest_range( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'belief', predicate('used_chest_range', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_used_chest_range', x};
end

function result = graph_stable_hr( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'belief', predicate('stable_hr', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_stable_hr', x};
end
function result = graph_stable_bf( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'belief', predicate('stable_bf', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_stable_bf', x};
end

% function result = graph_bel_phi( model, trace, parameters, t )
%   x   = l2.getall(trace, t+1, 'belief', predicate('phase_shift', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_bel_phi', x};
% end

% Support model

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
