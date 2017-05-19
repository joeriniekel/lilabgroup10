function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

%todo:
% dt 'vastzetten' op e.g. dt = 0.19

% adaptiemodel
% max en min zijn parameters, maar worden dus aangepast

% weerstand kan veranderen
% adaptiemodel kijkt naar de min en max chest_c en past aan de hand daarvan de formule
% voor de chest_c current naar chest_c in cm aan
% dus max chest_c is misschien 80cm
% max = 80cm; min = 70cm;
% chest_range = 10cm
% spanningsverschil wordt geobserveerd
% standaardformule: chest_c = (voltage^2 - a) * b
  % range voltage = max - min
  % neemt de weerstand kwadratisch toe?
  % a = minimum voltage (geen rek)
  % b = 'relatieve rek'
  %   b = range chest_c in cm / range voltage


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
  dir = trace(t).support.arg{1};
  % dir = '4 none';
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
  % sitfac        = trace(t).sitfac;

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
  % in bpm, between 0 and 180
  % output = hr * dt
  ps      = trace(t).ps.arg{1}; % t+1 doesn't work with this syntax + scenario values
  var     = trace(t+1).hr_var.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  a       = model.parameters.default.anxiety_hr;
  bhr     = model.parameters.default.bhr;
  lhr     = model.parameters.default.lhr;

  hr = (bhr + ps) + (a * anxiety) + var;
  if hr < lhr, hr = lhr; disp('lhr reached!'); end;

  global training_hr;
  if training_hr,    hr = 0;  end;
  result = {t+1, 'hr', hr};
end

function result = hr_pos( model, trace, parameters, t )
  %the 'input-current'
  global TRAINING;
  if TRAINING
    global TRAINING_HR
    pos = TRAINING_HR(t);
  else
    pos = 0;
  end
  result = {t+1, 'hr_pos', pos};
end


function result = breathing_f( model, trace, parameters, t )
  % between 0 and 4 (params.high_bf)
  hr_bpm  = trace(t+1).hr.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  h       = model.parameters.default.hr_breathing;
  % h2      = model.parameters.default.hr_breathing_exp;
  a2      = model.parameters.default.anxiety_bf;

  hr = hr_bpm / 60; % in s-1
  breathing_f = h * hr + (a2 * anxiety);

  global TRAINING;
  if TRAINING,    breathing_f = 0;  end; %bypass the domain model

  result = {t+1, 'breathing_f', breathing_f};
end

function result = used_chest_range( model, trace, parameters, t )
  % range = [0;1]
  % if higher than 1; subject is breathing moreintensely than healthy
  br_intensity = trace(t).breathing_intensity.arg{1};  %t+1 doesn't work for scenario predicates
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
  global CHEST_Y1 RT_CHEST1 %CHEST_Y2  RT_CHEST2
  CHEST_Y1(1) = [];
  CHEST_Y1(end+1) = curr_chest_c - 50;
  % CHEST_Y2(end+1) = curr_chest_c;
  refreshdata(RT_CHEST1);
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
  hr = trace(t+1).hr.arg{1};
  result = {t+1, 'observe', predicate('hr',hr)};
end
function result = obs_hr_pos( model, trace, parameters, t )
  pos = trace(t+1).hr_pos.arg{1};
  result = {t+1, 'observe', predicate('hr_pos',pos)};
end

function result = obs_chest_c( model, trace, parameters, t )
  chest_c = trace(t+1).chest_c.arg{1};
  result = {t+1, 'observe', predicate('chest_c',chest_c)};
end

function result = bel_chest_c( model, trace, parameters, t )
  chest_c = l2.getall(trace, t+1, 'observe', predicate('chest_c', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('chest_c',chest_c)};
end

function result = bel_starting_dir( model, trace, parameters, t )
  prev_chest_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  curr_chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  margin       = model.parameters.default.chest_c_margin;
  if curr_chest_c > prev_chest_c + margin
    chest_dir = '1 in';
  elseif curr_chest_c < prev_chest_c - margin
  	chest_dir = '3 out';
  else
  	chest_dir = '2 rest';
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
  prev_val = l2.getall(trace, a+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};

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
function result = bel_hr_pos( model, trace, parameters, t )
  hr_pos = trace(t+1).hr_pos.arg{1};
  result = {t+1, 'belief', predicate('hr_pos',hr_pos)};
end

%new
function result = bel_hr( model, trace, parameters, t )
  global training
  if ~training
    hr = l2.getall(trace, t+1, 'observe', predicate('hr', NaN)).arg{1}.arg{1};
  else
    % calculate the believed heart rate
    % this function functions the same as bel_breathing_f
    % search backwards; search until enough cycles are found
    n_cycles = model.parameters.default.n_hr_cycles;
    dt  = model.parameters.default.dt;
    max = 1/dt * model.parameters.default.max_interval_t; % same as bel bf
    floor = 1.0;

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
        if val > floor % make the value binary
          val = 1;
        else
          val = 0;
        end
        % modes:
        % if last (newest) observation = 0; mode = 0
        % if last (newest) observation > floor; mode = 1
        if isempty(mode)
          if val ~= prev_val
            mode = val; % 1 or 2
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

% bel breathing intensity



% %new
% function result = bel_relative_c( model, trace, parameters, t )
%   %the previous value of chest_c, relative to the used chest_range
%   chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
%   min          = model.parameters.default.min_chest_c; % believed min
%   max          = model.parameters.default.max_chest_c; % believed max
%   %todo param adaption gebruiken voor min en max?
%   range = max - min;
%   A = range / 2;
%   avg_chest_c = min + A;            % 70
%
%   relative_c = (chest_c - avg_chest_c) / A;
%     % value is in range [-1,1], except when the range changes between different points in time
%   if      relative_c >  1, relative_c =  1;
%   elseif  relative_c < -1, relative_c = -1;  end;
%
%   %ffor param adaption...
%   if chest_c > max - 0.5
%     relative_c = 1;
%   elseif chest_c < min + 0.5
%     relative_c = -1;
%   end
%
%   result = {t+1, 'belief', predicate('relative_c',relative_c)};
% end
% %new

% %new
% function result = bel_phase_shift( model, trace, parameters, t )
%   %phase_shift of chest_c-cycle on t-1,
%   % value will be used in (t-1) to calculate chest_c
%   % to calculate the phase, we need the new used_chest_range
%   starting_dir = l2.getall(trace, t+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
%   relative_c   = l2.getall(trace, t+1, 'belief', predicate('relative_c', NaN)).arg{1}.arg{1};
%   f            = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%   range = 1;  %assume that range is 1. if not, change min and max params?
%
%   % Formula for a basic sin wave: x(t) = A * sin(2*pi* f * t) + d;
%   % However, f can change during a cycle
%   % Adapt the formula so that f can be changed:
%   %   Let t = 1 and use +phi to shift the wave (with +phi in range[0,2pi])
%   %   phi has to be calculate manually because the range is also variable
%   %   (phi = asin( x(t)/prev_A) -2*pi*f*t) can not be used
%
%   % calculate the phase shift (phi)
%   %  2pi = full shift
%   %  1pi = half shift (inverted direction)
%   % .5pi = heighest point
%
%   % val = relative_c * 0.5 * pi;
%   % omschrijven: use asin to convert the [-1,1] to a 'sin-shape'
%   % want sin(t*0.5pi) is niet evenredig met t[0:1]
%   % val = asin(relative_c)/pi*2    * 0.5 * pi
%   val = asin(relative_c);
%
%   if strcmp(starting_dir, '1 in')
%     % phi = relative_c * 0.5 * pi;
%     phi = val;
%   else
%     %   % prev_starting_dir = '3 out'
%     % phi = relative_c * -0.5*pi + pi;
%     phi = val  * -1 + pi;
%     % eigenlijk:
%     % if relative_c > 0       % make pos relative_c negative
%     %   phi = val * -1 + pi;
%     % else                    % make negative relative_c positive
%     %   phi = val * -1 + pi;
%     % end
%   end
%   result = {t+1, 'belief', predicate('phase_shift',phi)};
% end
% %new
function result = bel_d_chest_c( model, trace, parameters, t )
  prev_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  c      = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('d_chest_c', c - prev_c)};
end
%new
function result = bel_avg_d_chest_c( model, trace, parameters, t )
  interval = 15;
  count = 0;
  v = [];
  if interval+1 > t, interval = t-1;end;
  for i=t:-1:t - interval
    count = count + 1;
    val = l2.getall(trace, i, 'belief', predicate('d_chest_c', NaN)).arg{1}.arg{1};
    v(count) = val;
  end
  % avg = mean(v);
  avg = 2;
  result = {t+1, 'belief', predicate('avg_d_chest_c', avg)};
end


% %new
function result = bel_d_hr( model, trace, parameters, t )
  prev_d  = l2.getall(trace, t, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
  prev_hr = l2.getall(trace, t, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  hr      = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  % laat verschil van t-1 ook meetellen: een stijging kan over meerdere tijdstappen plaatsvinden
  decay   = model.parameters.default.d_hr_decay;

  d = hr - prev_hr + prev_d * decay;
  result = {t+1, 'belief', predicate('d_hr', d)};
end

function result = bel_d_bf( model, trace, parameters, t )
  prev_d  = l2.getall(trace, t, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
  prev_bf = l2.getall(trace, t, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  bf      = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  decay   = model.parameters.default.d_hr_decay;

  d = bf - prev_bf + prev_d * decay;
  result = {t+1, 'belief', predicate('d_bf', d)};
end

%expected bf = hr * bf

function result = bel_anxiety( model, trace, parameters, t )
  %todo: niet naar d_bf kijken, maar naar expected_bf (from hr) vs belief_bf
  %new

  prev_anxiety = l2.getall(trace, t, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  bf    = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  hr    = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  % d_bf  = l2.getall(trace, t+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
  d_hr  = l2.getall(trace, t+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
  d_bf = 0;
  decay    = model.parameters.default.anxiety_decay;
  floor_hr = model.parameters.default.floor_hr;
  floor_bf = model.parameters.default.floor_bf;

  h       = model.parameters.default.hr_breathing;
  margin = 0;

  if onreglmatige d_intensity
  if onregelmatige d_bf
    anxiety = high
  elseif breathing intensity < hr * x
    anxiety = fairly high
  else
    anxiety = unknown
  end

  % d_bf + d_hr - komt door physical activity
  if d_bf + d_hr
    anxiety = unknown
  elseif d_bf && ~d_hr
    anxiety = high
  end





  if bf > h * hr - margin
    anxiety = 100 * (bf - h * hr) + prev_anxiety * decay;
    % anxiety = 10 * d_bf + prev_anxiety * decay;
  elseif d_hr > floor_hr
    % if hr is raised (by the physical state) anxiety cannot determined
    % thus it will stay at the same level
    anxiety = prev_anxiety;
    % disp('c1')
  elseif d_bf > floor_bf
    % disp('c2')
    anxiety = 10 * d_bf + prev_anxiety * decay;
  else
    % disp('c3')
    anxiety = prev_anxiety * decay;
  end
  % prev_anxiety
  % anxiety
  % anxiety = 0;
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
  floor_a = 0;
  if anxiety > floor_a
    assessment = true;
  else
    assessment = false;
  end
  result = {t+1, 'assessment', assessment};
end

% function result = bel_breathing_acc( model, trace, parameters, t )
%   % acceleration of the breathing f
%   % breathing f starts with 0, therefore acc = 0
%   breathing_f = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%   interval = model.parameters.default.acc_interval;
%   margin = model.parameters.default.acc_margin;
%
%   if interval > t-1, interval = t-1; end;
%   end_t = t - interval;
%   v = [];
%
%   for a=t:-1:end_t
%     val = l2.getall(trace, a, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%     count = t-a;              %teller: aantal loops
%     v(count+1) = val;
%     a = a-1;
%   end
%   %note that v is reversed, all the new (older t) values were appended to the vector
%
%   avg = sum(v) / length(v);
%   if length(v) == 0 disp('acc - div by 0'); end;
%   if breathing_f > avg + margin
%     breathing_acc = 1;
%   elseif breathing_f < avg - margin
%     breathing_acc = -1;
%   else
%     breathing_acc = 0;
%   end
%
%   result = {t+1, 'belief', predicate('breathing_acc',breathing_acc)};
% end
%
% function result = graph_bel_breathing_acc( model, trace, parameters, t )
%   acc = l2.getall(trace, t+1, 'belief', predicate('breathing_acc', NaN)).arg{1}.arg{1};
%   result = {t+1, 'graph_bel_breathing_acc', acc};
% end
%
% function result = bel_breathing_pattern( model, trace, parameters, t )
%
%   interval = model.parameters.default.pattern_interval;
%
%   if interval > t-1, interval = t-1; end;
%   end_t = t - interval;
%   v = [];
%
%   for a=t:-1:end_t
%     val = l2.getall(trace, a, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
%     count = t-a;              %teller: aantal loops
%     v(count+1) = val;
%     a = a-1;
%   end
%   %note that v is reversed, all the new (older t) values were appended to the vector
%
%   increases = v(v==1);
%   decreases = v(v==-1);
%   if ~isempty(increases) && ~isempty(increases)
%     pattern = 'irregular';
%   else
%     pattern = 'regular';
%   end
%
%   result = {t+1, 'belief', predicate('breathing_pattern',{pattern})};
% end
%
% function result = graph_bel_breathing_pattern( model, trace, parameters, t )
%   p = l2.getall(trace, t+1, 'belief', predicate('breathing_pattern', {NaN})).arg{1}.arg{1};
%   result = {t+1, 'graph_bel_breathing_pattern', {p}};
% end


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

% original_hr -> desire breathing_f
% belief phase shift + desire breathing_f -> desire starting direction

function result = des_bf( model, trace, parameters, t )
  % desired breathing_f without influence from anxiety.
  hr_bpm = l2.getall(trace, t+1, 'belief', predicate('original_hr', NaN)).arg{1}.arg{1};
  h       = model.parameters.default.hr_breathing;
  hr = hr_bpm / 60; % in s-1
  breathing_f = h * hr;

  %testing:   this should determine the final bf of the user
  % breathing_f = 0.4;
  result = {t+1, 'desire', predicate('breathing_f', breathing_f)};
end

% function result = des_starting_dir( model, trace, parameters, t )
%   % calculate the direction of the breathing cycle on a point t
%   % this is done by comparing x(t) and x(t+dt)
%   % to be precise, dt has to be as small as possible (lim->0) thus a new value for dt will be used.
%   % (t-dt) and (t+dt) is more realistic, but this can prevent the cycle from completing in certain instance.
%   t1 = 1;
%   t2 = 1.0001;
%   f   = l2.getall(trace, t, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
%   phi = l2.getall(trace, t, 'belief', predicate('phase_shift', NaN)).arg{1}.arg{1};
%   dt  = model.parameters.default.dt;
%   % x1 = A * sin(2*pi* f * t1 + phi) + avg_chest_c;
%   % x2 = A * sin(2*pi* f * t2 + phi) + avg_chest_c;
%   % A, avg_chest_c can be disregarded
%   x1 = sin(2*pi*f*dt* t1 + phi);
%   x2 = sin(2*pi*f*dt* t2 + phi);
%
%   if x1 < x2
%     curr_dir2 = '1 in';
%   else
%     curr_dir2 = '3 out';
%   end
%
%   result = {t+1, 'desire', predicate('starting_dir', curr_dir2)};
% end




% alt:
% if breathing_f < des bf
%   change dir
% if breathing_f > des bf
%   advice rest


% oude des_starting_dir werkt niet:
%   support is afhankelijk van chest c, maar stuurt de breathing intensity
%   (snelheid van het ademen zelf) niet aan.
%
%

%when to start cycle time? todo

function result = cycle_time( model, trace, parameters, t )
  % we can just give support about breathing in or out (quantitative)
  % not exactly how deep (fast) to breathe (qualitative)

  % when user should breathe in:  cycle_time = positive
  % when user should breathe out: cycle_time = negative
  % cycle_time is never 0

  %if assessment(t) = false and assessment(t+1) = true
  %   cycle time = 'relative_c'/100 * max_cycle_time
  prev_ct = trace(t).cycle_time.arg{1};
  bf      = l2.getall(trace, t+1, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
  dt  = model.parameters.default.dt;
  max_cycle_time = 1/bf * 1/dt; % 1/0.5 * 1/0.1 = 20

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
  result = {t+1, 'cycle_time', cycle_time};
end

function result = des_starting_dir( model, trace, parameters, t )
  cycle_time = trace(t+1).cycle_time.arg{1};
  chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  min          = model.parameters.default.min_chest_c; % believed min
  max          = model.parameters.default.max_chest_c; % believed max
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
% function result = des_starting_dir( model, trace, parameters, t )
%   % prevent giving a false desired starting dir when phi = 0.49pi
%   % compute average breathing direction in the next interval
%   f   = l2.getall(trace, t, 'desire', predicate('breathing_f', NaN)).arg{1}.arg{1};
%   phi = l2.getall(trace, t, 'belief', predicate('phase_shift', NaN)).arg{1}.arg{1};
%   dt = model.parameters.default.dt;
%
%   margin = asin(1);  %90%
%   % alsublieft, maar 1 enkel if-statement
%   % e.g. if phi in between [margin, pi - margin]; 'rest'
%   if (margin < phi && phi < pi - margin) || (pi + margin < phi && phi < 2*pi - margin)
%     avg_dir = '2 rest'
%     margin
%     phi
%   else
%     %params
%     start = 0.1;
%     interval = 0.3;
%     dirs = [];
%     count = 0;
%     for t1=start:interval:1
%       x1 = sin(2*pi*f*dt* t1 + phi);
%       x2 = sin(2*pi*f*dt* t1 + 0.001 + phi);
%       count = count+1;
%       if x1>=x2
%         dirs(count) = 1;
%       else
%         dirs(count) = -1;
%       end
%     end
%
%     if sum(dirs)/count > 0
%       avg_dir = '1 in';
%     else
%       avg_dir = '3 out';
%     end
%   end
%   result = {t+1, 'desire', predicate('starting_dir', {avg_dir})};
% end


function result = support( model, trace, parameters, t )
  starting_dir = l2.getall(trace, t+1, 'desire', predicate('starting_dir', NaN)).arg{1}.arg{1};
  assessment = trace(t+1).assessment.arg{1};
  % assessment = true;

  if assessment
    starting_dir = starting_dir;
  else
    starting_dir = '4 none';
  end
  % if t == 1, starting_dir =  '4 none';end;
  % starting_dir =  '4 none';
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
function result = graph_d_hr( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_d_hr', x};
end
function result = graph_d_bf( model, trace, parameters, t )
  x = l2.getall(trace, t+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_d_bf', x};
end

% Starting direction, physical state, original_hr, anxiety

function result = graph_bel_starting_dir( model, trace, parameters, t )
  chest_pos = l2.getall(trace, t+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
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
