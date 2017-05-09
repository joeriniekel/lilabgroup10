function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

% ---------------------------------------------------------------
%
%   DOMAIN
%
% ---------------------------------------------------------------

function result = anxiety_regulation( model, trace, parameters, t )
  breathing_f = trace(t).breathing_f.arg{1};
  disfac = model.parameters.default.disfac_on_regulation;
  b = model.parameters.default.breathing_anxiety;

  anxiety_regulation = disfac + b * (1 - breathing_f);

  result = {t+1, 'anxiety_regulation', anxiety_regulation};
end

function result = anxiety( model, trace, parameters, t )
  prev_anxiety  = trace(t).anxiety.arg{1};
  sitfac        = trace(t).sitfac.arg{1}; %gives errors (t+1 of t)
  regulation    = trace(t+1).anxiety_regulation.arg{1};
  decay         = model.parameters.default.anxiety_decay;
  disfac        = model.parameters.default.disfac;
  s             = model.parameters.default.sitfac_anxiety;

  % breathing pattern
  % inconsistent

  new_anxiety = (s * sitfac + prev_anxiety * decay) * regulation;
  result = {t+1, 'anxiety', new_anxiety};
end

function result = hr( model, trace, parameters, t )
  dt      = model.parameters.default.dt;
  ps      = trace(t).physical_state.arg{1}; % t+1 doesn't work with this syntax + scenario values
  anxiety = trace(t+1).anxiety.arg{1};
  a       = model.parameters.default.anxiety_hr;
  bhr     = dt * model.parameters.default.bhr; %todo kan dit?
  lhr     = dt * model.parameters.default.lhr;

  hr_new = (bhr * ps) + (a * anxiety) + 0.01 * rand;
  if hr_new < lhr, hr_new = lhr; disp('lhr reached!'); end;

  result = {t+1, 'hr', hr_new};
end

function result = breathing_f( model, trace, parameters, t )
  hr  = trace(t+1).hr.arg{1};
  h   = model.parameters.default.hr_breathing;
  h2  = model.parameters.default.hr_breathing_exp;

  breathing_f_in_bpm = h * hr ^ h2;
  breathing_f = breathing_f_in_bpm / 60;% + 0.001*rand;
  result = {t+1, 'breathing_f', breathing_f};
end

function result = used_chest_range( model, trace, parameters, t )
  br_intensity = trace(t).breathing_intensity.arg{1};  %t+1 doesn't work for scenario predicates
  max = model.parameters.default.max_chest_range;

  used_chest_range = max * br_intensity;
  result = {t+1, 'used_chest_range', used_chest_range};
end

function result = prev_relative_c( model, trace, parameters, t )
  %the value of chest_c, relative to the preferred/used chest_range
  prev_chest_c      = trace(t).chest_c.arg{1};           % 60-80
  range             = trace(t+1).used_chest_range.arg{1}; % 19
  min               = model.parameters.default.min_chest_c;
  max               = model.parameters.default.max_chest_c;
  avg_chest_c       = min + ((max - min) / 2);            % 70
  A = range / 2;

  relative_c = (prev_chest_c - avg_chest_c) / A;
    % value is in range [-1,1], except when the range changes during simulation
  if      relative_c >  1, relative_c =  1;
  elseif  relative_c < -1, relative_c = -1;  end;

  result = {t+1, 'prev_relative_c', relative_c};
end

function result = prev_phase_shift( model, trace, parameters, t )
  %phase_shift of chest_c-cycle on t-1,
  % value will be used in (t-1) to calculate chest_c
  % to calculate the phase, we need the new used_chest_range

  prev_starting_dir = trace(t).starting_dir.arg{1};
  relative_c        = trace(t+1).prev_relative_c.arg{1};
  range             = trace(t+1).used_chest_range.arg{1};
  f                 = trace(t+1).breathing_f.arg{1};

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

  %%% use asin to convert the [-1,1] to a 'sin-shape'   (0.25 and 0.75 are mostly affected)
  %%% val = asin(val)?
  if strcmp(prev_starting_dir, '1 in')
    phi = relative_c * 0.5 * pi;
  else
    %   % prev_starting_dir = '3 out'
    phi = relative_c * -0.5*pi + pi;
    % eigenlijk:
    % if relative_c > 0       % make pos relative_c negative
    %   phi = val * -1 + pi;
    % else                    % make negative relative_c positive
    %   phi = val * -1 + pi;
    % end
  end
  result = {t+1, 'prev_phase_shift', phi}; % + 0.25*pi
end

function result = chest_c( model, trace, parameters, t )
  prev_chest_c      = trace(t).chest_c.arg{1};

  f     = trace(t+1).breathing_f.arg{1};
  range = trace(t+1).used_chest_range.arg{1};
  phi   = trace(t+1).prev_phase_shift.arg{1};

  dt    = model.parameters.default.dt;
  min   = model.parameters.default.min_chest_c;
  max   = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2); % 80
  % f = dt * f;
  A = range / 2;
  % phi = phi * -1;
  curr_chest_c = A * sin(2*pi* f * dt + phi) + avg_chest_c;

  %oude fallback
  margin = 0.2;
  add = 0;
  if curr_chest_c == prev_chest_c + margin
    disp('curr_chest_c == prev_chest_c')
    curr_chest_c = curr_chest_c + add;
  elseif curr_chest_c == prev_chest_c - margin
    disp('curr_chest_c == prev_chest_c')
    curr_chest_c = curr_chest_c + add;
  end
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
  phi   = trace(t+1).prev_phase_shift.arg{1};
  f     = trace(t+1).breathing_f.arg{1};
  dt    = model.parameters.default.dt;
  % x1 = A * sin(2*pi* f * t1 + phi) + avg_chest_c;
  % x2 = A * sin(2*pi* f * t2 + phi) + avg_chest_c;
  % A, f, avg_chest_c can be disregarded
  x1 = sin(2*pi*f*dt* t1 + phi);
  x2 = sin(2*pi*f*dt* t2 + phi);

  if x1 < x2, curr_dir2 = '1 in';
  else
    curr_dir2 = '3 out';
  end

  result = {t+1, 'starting_dir', {curr_dir2}};
end

















% ---------------------------------------------------------------
%
%   ANALYSIS
%
% ---------------------------------------------------------------



function result = obs_chest_c( model, trace, parameters, t )
  chest_c = trace(t).chest_c.arg{1};
  % chest_c = l2.getall(trace, t, 'chest_c', {NaN}).arg{1};
  result = {t+1, 'observe', predicate('chest_c',chest_c)};
end

function result = bel_chest_c( model, trace, parameters, t )
  chest_c = l2.getall(trace, t+1, 'observe', predicate('chest_c', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('chest_c',chest_c)};
end

function result = bel_chest_pos( model, trace, parameters, t )
  curr_chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  prev_chest_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  margin       = model.parameters.default.chest_c_margin;

  if curr_chest_c > prev_chest_c + margin
    chest_pos = '1 in';
  elseif curr_chest_c < prev_chest_c - margin
  	chest_pos = '3 out';
  else
  	chest_pos = '2 rest';
  end
  result = {t+1, 'belief', predicate('chest_pos',{chest_pos})};
end

function result = graph_bel_chest_pos( model, trace, parameters, t )
  chest_pos = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_chest_pos', {chest_pos}};
end


function result = bel_chest_change( model, trace, parameters, t )
  %check whether a breathing cycle has started
  prev_pos = l2.getall(trace, t, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  curr_pos = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  %param om in- of uitademen te meten
  mode = model.parameters.default.chest_change_mode;

  %check if there is a change,
  % when there is, check if the breathing cycle has restarted
  if strcmp(curr_pos,prev_pos) %no change
    chest_change = false;
  elseif strcmp(curr_pos,'1 in') && mode == 1
    chest_change = true;
  elseif strcmp(curr_pos,'3 out') && mode == 3
    chest_change = true;
  else
    chest_change = false;
  end
  result = {t+1, 'belief', predicate('chest_change',chest_change)};
end

function result = graph_bel_chest_change( model, trace, parameters, t )
  c = l2.getall(trace, t+1, 'belief', predicate('chest_change', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_chest_change', c};
end

function result = bel_breathing_f( model, trace, parameters, t )
  %calculate believed breathing frequency
  %search backwards; search until enough breathing cycles are found
  n_cycles = 3;
  max = 10000;  %afhankelijk van params maken --> dt...

  n_cycles = model.parameters.default.n_breathing_cycles;
  max = model.parameters.default.max_interval_t;
  dt = model.parameters.default.dt;

  max = max / dt;
  v = [];
  a = t;

  if t < n_cycles   % a minimum of n*2 time steps are required
    breathing_f = 0;
    disp('Waiting - collecting data')
  else
    if t==n_cycles, disp('  Go  -->'); end;
    while sum(v) <= n_cycles
      val = l2.getall(trace, a+1, 'belief', predicate('chest_change', NaN)).arg{1}.arg{1};
      count = t-a;              %teller: aantal loops
      v(count+1) = val;         %matlab starts counting at 1 instead of 0
      a = a-1;
      if a <= 1,    break; end; %break at t1
      if count>max, break; end; %break when max timesteps has been searched
    end
    %note that v is reversed, all the new (older t) values were appended to the vector

    %calulate interval size + number of breathing cycles
    t_last_start = t - find(v==1,1); %can be [] (empty)
    if isempty(t_last_start)          %no cycles found: return 0
      breathing_f = 0;
      % disp('no cycles found')
    else
      t_start = a;
      t_end = t_last_start - 1;
      if t_end<1, t_end=1; disp(t); disp('t_end < 1'); end;
      if t_end<t_start, disp('kanniet - interval onjuist');  end;

      interval = t_end - t_start;
      n_cycles = sum(v);    %can be less than the original
      breathing_f = n_cycles / interval;
    end
  end

  result = {t+1, 'belief', predicate('breathing_f', breathing_f)};
end

function result = graph_bel_breathing_f( model, trace, parameters, t )
  b = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_breathing_f', b};
end
function result = graph_breathing_f_error( model, trace, parameters, t )
  b = l2.getall(trace, t+1, 'breathing_f', NaN).arg{1};
  c = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_breathing_f_error', c-b};
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
