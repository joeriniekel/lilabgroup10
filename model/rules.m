function [ fncs ] = rules()    % DO NOT EDIt
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

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
  % [0;100]
  sitfac        = trace(t).sitfac.arg{1}; %gives errors (t+1 of t)
  prev_anxiety  = trace(t).anxiety.arg{1};
  regulation    = trace(t+1).anxiety_regulation.arg{1};
  decay         = model.parameters.default.anxiety_decay;
  disfac        = model.parameters.default.disfac;
  s             = model.parameters.default.sitfac_anxiety;
  % Inf waardes geven errors... wachten op bugfix van linford
  % sitfac        = trace(t).sitfac;
  % t
  % sitfac
  % sitfac = 0;
  new_anxiety = (s * sitfac + prev_anxiety * decay) * (1 - regulation);
  result = {t+1, 'anxiety', new_anxiety};
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

  hr_new = (bhr + ps) + (a * anxiety) + var;
  if hr_new < lhr, hr_new = lhr; disp('lhr reached!'); end;
  result = {t+1, 'hr', hr_new};
end

function result = breathing_f( model, trace, parameters, t )
  % between 0 and 4 (params.high_bf)
  hr_bpm  = trace(t+1).hr.arg{1};
  anxiety = trace(t+1).anxiety.arg{1};
  h       = model.parameters.default.hr_breathing;
  h2      = model.parameters.default.hr_breathing_exp;
  a2      = model.parameters.default.anxiety_bf;

  hr = hr_bpm / 60; % in s-1
  breathing_f = h * hr + (a2 * anxiety);
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

function result = prev_relative_c( model, trace, parameters, t )
  %the previous value of chest_c, relative to the preferred/used chest_range
  prev_chest_c = trace(t).chest_c.arg{1};            % 60-80
  range        = trace(t+1).used_chest_range.arg{1}; % 19
  min          = model.parameters.default.min_chest_c;
  max          = model.parameters.default.max_chest_c;
  avg_chest_c  = min + ((max - min) / 2);            % 70
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

  % val = relative_c * 0.5 * pi;
  % omschrijven: use asin to convert the [-1,1] to a 'sin-shape'
  % want sin(t*0.5pi) is niet evenredig met t[0:1]
  % val = asin(relative_c)/pi*2    * 0.5 * pi
  val = asin(relative_c); %deze regel voorkomt vastlopen van de grafiek op relative_c = 1 of -1

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
  result = {t+1, 'prev_phase_shift', phi}; % + 0.25*pi
end

function result = chest_c( model, trace, parameters, t )
  prev_chest_c      = trace(t).chest_c.arg{1};

  f     = trace(t+1).breathing_f.arg{1};
  range = trace(t+1).used_chest_range.arg{1};
  phi   = trace(t+1).prev_phase_shift.arg{1};

  dt          = model.parameters.default.dt;
  min         = model.parameters.default.min_chest_c;
  max         = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2); % 80
  % f = dt * f; % breathing_f is already inheritably adapted to dt,
                % but if dt is changed during simulation, it could produce weird behaviour.
  A = range / 2;
  % y(t) = A sin(2 pi f t dt + phi)
  curr_chest_c = A * sin(2*pi* f * dt + phi) + avg_chest_c;
  %oude fallback
  % margin = 0.2;
  % add = 0;
  % if curr_chest_c >= prev_chest_c + margin  %(and direction is 'in')
  %   disp('curr_chest_c == prev_chest_c')
  %   curr_chest_c = curr_chest_c - add;
  % elseif curr_chest_c <= prev_chest_c - margin
  %   disp('curr_chest_c == prev_chest_c')
  %   curr_chest_c = curr_chest_c + add;
  % end
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
  %dt...
  result = {t+1, 'observe', predicate('hr',hr)};
end

function result = obs_chest_c( model, trace, parameters, t )
  chest_c = trace(t+1).chest_c.arg{1};
  result = {t+1, 'observe', predicate('chest_c',chest_c)};
end

function result = bel_hr( model, trace, parameters, t )
  hr = l2.getall(trace, t+1, 'observe', predicate('hr', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('hr',hr)};
end

function result = bel_chest_c( model, trace, parameters, t )
  chest_c = l2.getall(trace, t+1, 'observe', predicate('chest_c', NaN)).arg{1}.arg{1};
  result = {t+1, 'belief', predicate('chest_c',chest_c)};
end

function result = bel_starting_dir( model, trace, parameters, t )
  curr_chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  prev_chest_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
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
  mode     = model.parameters.default.breathing_dir_mode;
  max = max / dt;   % t=1,60x     t=0.1,600x
  v = [];
  a = t;
  prev_val = '';

  if t <= n_cycles*2   % a minimum of n*2+1 time steps are required (in+out+...+next_in)
    breathing_f = 0;
    disp('Waiting - collecting data')
    if t==n_cycles*2, disp('  Go  -->'); end; %last wait
  else
    while sum(v) <= n_cycles + 1      % the last interval will be cut off
      val = l2.getall(trace, a+1, 'belief', predicate('starting_dir', NaN)).arg{1}.arg{1};
      count = t-a;              %teller: aantal loops
      if (strcmp(val, '1 in') && mode==1) || (strcmp(val, '3 out') && mode==3)
        if ~strcmp(val,prev_val)  %check if the cycle has started
          v(count+1) = 1;         %matlab starts counting at 1 instead of 0
        else
          v(count+1) = 0;
        end
      else
        v(count+1) = 0;
      end
      prev_val = val;
      a = a-1;
      if a <= 1,      break; end; %break at t1
      if count>max,   break; end; %break when max timesteps has been searched
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
      if t_end<1 && t>10, t_end=1; disp(t); disp('t_end < 1'); end;
      if t_end<t_start, disp('kanniet - interval onjuist. t =');disp(t);  end;

      interval = t_end - t_start;
      n = sum(v);    %can be less than the original
      breathing_f = n / interval * dt;
    end
  end
  result = {t+1, 'belief', predicate('breathing_f', breathing_f)};
end

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

function result = bel_anxiety( model, trace, parameters, t )
  prev_anxiety = l2.getall(trace, t, 'belief', predicate('anxiety', NaN)).arg{1}.arg{1};
  bf = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  hr = l2.getall(trace, t+1, 'belief', predicate('hr', NaN)).arg{1}.arg{1};
  d_bf = l2.getall(trace, t+1, 'belief', predicate('d_bf', NaN)).arg{1}.arg{1};
  d_hr = l2.getall(trace, t+1, 'belief', predicate('d_hr', NaN)).arg{1}.arg{1};
  decay     = model.parameters.default.anxiety_decay;
  floor_hr  = 15;
  floor_bf  = 0.05;
  if d_hr > floor_hr
    % if hr is raised (by the physical state) anxiety cannot determined
    % thus it will stay at the same level
    anxiety = prev_anxiety;
    disp('c1')
  elseif d_bf > floor_bf
    disp('c2')
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

function result = bel_prev_ps( model, trace, parameters, t )
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
  floor_a = 10;
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






% Breathing_f and heart rate

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
