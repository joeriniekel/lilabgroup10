function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end


% --------------------------------------------------------
%
%   DOMAIN
%
% --------------------------------------------------------


function result = anxiety_regulation( model, trace, parameters, t )
  breathing_f = trace(t).breathing_f.arg{1};
  disfac = model.parameters.default.disfac_on_regulation;
  b = model.parameters.default.breathing_anxiety;

  anxiety_regulation = disfac + b * (1 - breathing_f);

  result = {t+1, 'anxiety_regulation', anxiety_regulation};
end

% oude anxiety formule
% function result = anxiety( model, trace, parameters, t )
%   % sitfac2 = trace(t).sitfac.arg{1}; %gives errors
%   for sf = l2.getall(trace, t, 'sitfac', {NaN}) %t+1? todo
%     sitfac = sf.arg{1};
%
%     prev_anxiety = trace(t).anxiety.arg{1};
%     %sitfac = trace(t).sitfac.arg{1}; %gives errors
%     breathing_f = trace(t).breathing_f.arg{1};
%
%     decay = model.parameters.default.anxiety_decay;
%     disfac = model.parameters.default.disfac;
%     s = model.parameters.default.sitfac_anxiety;
%     b = model.parameters.default.breathing_anxiety;
%     %dispositional factors affect the intensity of anxiety provoking
%     %situations
%     new_anxiety = prev_anxiety * decay + ...
%       disfac * (s * sitfac + b * breathing_f);
%
%     result = {t+1, 'anxiety', new_anxiety};
%   end
% end
function result = anxiety( model, trace, parameters, t )
  prev_anxiety = trace(t).anxiety.arg{1};
  sitfac = trace(t).sitfac.arg{1}; %gives errors (t+1 of t)
  regulation = trace(t+1).anxiety_regulation.arg{1};

  decay = model.parameters.default.anxiety_decay;
  disfac = model.parameters.default.disfac;
  s = model.parameters.default.sitfac_anxiety;
    %dispositional factors affect the intensity of anxiety provoking
    %situations
    % new_anxiety = prev_anxiety * decay + ...
    %   disfac * (s * sitfac + b * breathing_f);
  new_anxiety = (s * sitfac + prev_anxiety * decay) * regulation;

  result = {t+1, 'anxiety', new_anxiety};
end

function result = hr( model, trace, parameters, t )

  dt = model.parameters.default.dt;
  anxiety = trace(t+1).anxiety.arg{1};
  bhr = dt * model.parameters.default.bhr;
  ps = trace(t).physical_state.arg{1}; % t+1 doesn't work with this syntax + scenario values
  a = model.parameters.default.anxiety_hr;

  hr_new = (bhr * ps) + (a * anxiety);  % * rand

  result = {t+1, 'hr', hr_new};
end

function result = breathing_f( model, trace, parameters, t )

  hr = trace(t+1).hr.arg{1};
  h = model.parameters.default.hr_breathing;
  h2 = model.parameters.default.hr_breathing_exp;

  breathing_f_in_bpm = h * hr ^ h2;
  breathing_f = breathing_f_in_bpm / 60;% + 0.001*rand;

  breathing_f = 0.21; %todo

  result = {t+1, 'breathing_f', breathing_f};
end





% breathing_intensity
% the relation of this concept is unknown. it should look like this:
%
% function result = breathing_intensity( model, trace, parameters, t )
%   % hoe sneller iemand ademt, hoe minder diep elke ademhaling is
%   % bij snel ademen, ademt iemand maar voor 50% in.
%   breathing_f = trace(t+1).breathing_f.arg{1};
%   max_breathing_f = 3;
%   % min_breathing_f = 0;
%   breathing_intensity = breathing_f / max_breathing_f;
%
%   result = {t+1, 'breathing_intensity', breathing_intensity};
% end

function result = used_chest_range( model, trace, parameters, t )
  breathing_intensity = trace(t).breathing_intensity.arg{1};  %t+1 doesn't work for scenario predicates
  max = model.parameters.default.max_chest_range;
  used_chest_range = max * breathing_intensity;
  result = {t+1, 'used_chest_range', used_chest_range};
end

%nutteloos concept
% function result = chest_speed( model, trace, parameters, t )
%   breathing_f = trace(t+1).breathing_f.arg{1};
%   used_chest_range = trace(t+1).used_chest_range.arg{1};
%   chest_speed = used_chest_range / breathing_f;
%   result = {t+1, 'chest_speed', chest_speed};
% end

% function result = chest_pos3( model, trace, parameters, t )
% %  prev_chest_c = trace(t).chest_c3.arg{1};
%   if t>2
%     curr_chest_c = trace(t+1).chest_c3.arg{1};
%   else
%     curr_chest_c = prev_chest_c;
%   end
%   if prev_chest_c < curr_chest_c
%     %agent was breathing in
%     %he will continue to do so if max range is not reached in dt
%     phi = 0;
%   else
%     %agent was breathing out
%     phi = 0.5*pi;   %pi is exactly half a cycle
%   end
% % end



function result = chest_c3( model, trace, parameters, t )
  prev_chest_c = trace(t).chest_c3.arg{1};
  if t>2
    curr_chest_c = trace(t+1).chest_c3.arg{1};
  else
    curr_chest_c = prev_chest_c;
  end

  %t+2...
  range = trace(t+1).used_chest_range.arg{1};
  f = trace(t+1).breathing_f.arg{1};
  dt = model.parameters.default.dt;
  min = model.parameters.default.min_chest_c;
  max = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2);
  A = range / 2;

  relative_c = (curr_chest_c - avg_chest_c) / range;
    % [-1,1]
  if prev_chest_c < curr_chest_c
    %starting direcion = breathe in
    if relative_c >= 0
      phase = relative_c * 0.5*pi;
    else
      phase = relative_c * 0.5*pi; % (-)
    end
  else
    %starting direcion = breathe out
    if relative_c >= 0
      phase = relative_c * 0.5*pi * -1 + pi;
    else
      phase = relative_c * 0.5*pi * -1 + pi; % (-)
    end
  end

  if relative_c > 1, disp(relative_c); end;
  if relative_c == 0, disp(relative_c); end;
  % phase = relative_c * 0.5*pi + phi;
  % phase = relative_c * polarity * 0.5*pi + add_pi*pi;

  % standard formula
  % x(t) = A * sin(2*pi*f*t) + d
  % d = default level = avg_chest_c
  % A = amplitude = max available range + or - avg_chest_c
  % here; t = dt, not the actual t
  % x(dt) = A * sin(2*pi*f*dt + phase) + avg_chest_c
  % si2n = sin(2*pi*f*dt)

  % phase
  % aa = A * sin(2*pi*f*dt + phase)
  new_chest_c = A * sin(2*pi*f*dt + phase) + avg_chest_c;

  result = {t+2, 'chest_c3', new_chest_c};
end

function result = chest_pos_phi( model, trace, parameters, t )
  %determine the chest position/transition by calculating the phase
  % range = [0, 2pi]
  %first a new waveform will be recalculated with a different phase for a different f
  %this can be used to calculate the new x(t) and thus the new chest_c

  %note that phi itself will have a sawtooth waveshape instead of a sin
  dt = model.parameters.default.dt;
  sig = model.parameters.default.sig;   %significance
  lower_bound = model.parameters.default.lower_bound;
  %previous values
  prev_f = trace(t).breathing_f.arg{1};
  prev_t = t - 1;
  prev_phi = trace(t).chest_pos_phi.arg{1};
  prev_range = trace(t).used_chest_range.arg{1};
  prev_A = prev_range / 2;
  %current values
  f = trace(t+1).breathing_f.arg{1};
  range = trace(t+1).used_chest_range.arg{1};
  A = range / 2;

  % prediction for x(t+1),     (because t is filled in...)
  % x1 = x(prev_t) = prev_A * sin(2*pi * prev_f * prev_t + 0);
  x1 = prev_A * sin(2*pi * prev_f * prev_t * dt + prev_phi);
  x1 = round(x1,sig);
  if x1 < lower_bound, x1 = lower_bound; end;      %ignore tiny values
  %because f can change,
  % x2 = x(prev_t) = prev_A * sin(2*pi * f * prev_t + phi);
  % Thus: x1 = x2 = x(prev_t),  but x2 has a different breathing f
  % phi is unknown; subtract phi:
  % new formula - omschrijven:
  % x2 = x1 = prev_A * sin(2*pi*new_f*prev_t + new_phi);
  % 2*pi*new_f*prev_t + new_phi = asin( x(prev_t) / A) )
  % new_phi = asin( x(prev_t) / A) ) - 2*pi*new_f*prev_t
  phi = asin( x1 / prev_A) - 2*pi * f * prev_t * dt;
  phi = mod(phi,2*pi);
  % prev_f
  % f
  x3 = prev_A * sin(2*pi * f * prev_t * dt + phi);
  x4 = prev_A * sin(2*pi * prev_f * prev_t * dt + phi);
  % disp(x4 - x3)
  % d_phi = (f * dt) * 2*pi;
  % phi = phi + d_phi;
  % phi = 2*pi - phi;

  % x2 should have the same value as x1, even though the f is changed
  % x2 = x(prev_t) = prev_A * sin(2*pi * f * prev_t * dt + phi);
  %%%%    code om het te testen
  %%%% x2 = prev_A * sin(2*pi * f * prev_t * dt + phi);
  %%%% x2 = round(x2,sig);
  %%%% if x2 < lower_bound, x2 = lower_bound; end;
  % phi = 0;
  result = {t+1, 'chest_pos_phi', phi};
end

%alternatieve manier

% chest_c veroorzaakt nu schokken
% een betere modulatie zou zijn:
%   bepaal huidige transition:  in- of uitademen  (in dit geval fase)
%   bepaal breathing f
%   bepaal de snelheid waarmee de chest pos moet bewegen
%     (f/afstand)
%   bepaal aan de hand daarvan de nieuwe chest pos (e.g. oude pos + 2cm)
%
%   op deze manier is het modeleren van afwijkingen ook logischer



function result = chest_c( model, trace, parameters, t )
  dt = model.parameters.default.dt;
  sig = model.parameters.default.sig;   %significance
  lower_bound = model.parameters.default.lower_bound;
  phi = trace(t+1).chest_pos_phi.arg{1};  %t+1 ?

  f = trace(t+1).breathing_f.arg{1};
  range = trace(t+1).used_chest_range.arg{1};
  A = range / 2;
%---------------------------------------%
                                        %
  % twee opties                         %
  x2 = A * sin(2*pi* f * 1 * dt + phi); %
  % of                                  %
  x2 = A * sin(2*pi* f * t * dt + 0);   %
%---------------------------------------%
  x2 = round(x2,sig);
  if x2 < lower_bound, x2 = lower_bound; disp('lb'); disp(t); end;

  min = model.parameters.default.min_chest_c;
  max = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2);

  chest_c = avg_chest_c + x2;

  result = {t+1, 'chest_c', x2};
end

% oude chest_c
function result = chest_c2( model, trace, parameters, t )
    dt = model.parameters.default.dt;
    max_c = model.parameters.default.max_chest_c;
    min_c = model.parameters.default.min_chest_c;
    breathing_f = trace(t+1).breathing_f.arg{1};
    %breathing_f = n/s
    %n = breathing_f * s;

    n = breathing_f * t * dt;   %number of breathing cycles
    phase = mod(n,1);     %the phase of the cycle: a number between 0 and 1
    %in this case the breathing cycle always starts at t=0

    if phase < 0.5
        phi = phase;
    else
        phi = 1 - phase;
    end
    chest_c = max_c * 2 * phi;

    result = {t+1, 'chest_c2', chest_c};
end

% oude chest_pos
function result = chest_pos( model, trace, parameters, t )
  curr_chest_c = trace(t+1).chest_c.arg{1};
  prev_chest_c = trace(t).chest_c.arg{1};
  margin = model.parameters.default.margin;

  if curr_chest_c > prev_chest_c + margin
    chest_pos = '1 in';
  elseif curr_chest_c < prev_chest_c - margin
   	chest_pos = '3 out';
  else
   	chest_pos = '2 rest';
  end
  result = {t+2, 'chest_pos', {chest_pos}};
end


function result = prev_phase_shift( model, trace, parameters, t )
  %phase_shift of chest_c-cycle on t-1,
  % value will be used in t to calculate chest_c
  % to calculate the phase, we need the new used_chest_range
  t
  prev_chest_c      = trace(t).chest_c4.arg{1};
  prev_starting_dir = trace(t).starting_dir.arg{1};
  range             = trace(t+1).used_chest_range.arg{1};
  f                 = trace(t+1).breathing_f.arg{1};

  dt                = model.parameters.default.dt;
  min               = model.parameters.default.min_chest_c;
  max               = model.parameters.default.max_chest_c;
  avg_chest_c       = min + ((max - min) / 2);

  %formula for a basic sin wave:  x(t) = A * sin(2*pi* f * t) + d;
  % This: chest_c = A * sin(2*pi* f * t) + avg_chest_c;
  % However, f can change during a cycle
  % Adapt the formula so that f can be changed
  %   Let t = 1 and use +phi to shift the wave (with +phi in range[0,2pi])
  %   phi has to be calculate manually because the range is also variable
  %   (phi = asin( x(t)/prev_A) -2*pi*f*t) can not be used
  % range = range +10;
  relative_c = (prev_chest_c - avg_chest_c) / range
    % value is in range [-1,1], except when the range changes
  if      relative_c >  1, relative_c =  1;
  elseif  relative_c < -1, relative_c = -1;  end;
  % relative_c = asin(relative_c);

  % calculate the phase shift (phi)
  % 2pi = a full shift
  % 1pi = a half shift (inverted direction)

  % prev_starting_dir = '1 in';

  % use asin to convert the [-1,1] to a 'sin-shape'   (0.25 and 0.75 are mostly affected)
  % phi = asin(relative_c) * pi   %convert the value to pi-rad

  if strcmp(prev_starting_dir, '1 in')
    %convert to 2pi
    phi = relative_c * 0.5 * pi
  else
    %   % prev_starting_dir = '3 out'
    if relative_c > 0
      phi = relative_c * -1 * 0.5 * pi
    else
      phi = relative_c * -1 * 0.5 * pi
    end
  end

  % plot(t,sin(2*pi*f*t + 0.2*pi + pi))

  % if strcmp(prev_starting_dir, '1 in')
  %   phi = relative_c * 0.5*pi;
  %   phi2 = acos(relative_c);
  %   phi3 = asin( relative_c ) - 2*pi*f;
  %   % 20% = 0.2 corresponds to a phase shift of   +0.1pi
  %   %      -0.2   ~   -0.1pi or + (2 - 0.1)pi = 1.9pi
  %     %   sin(2*pi - (2 - 0.1)*pi) = sin(2*pi - 0.1*pi + pi)
  % else
  %   % prev_starting_dir = '3 out'
  %   phi = relative_c * -0.5*pi + pi;
  %   phi2 = acos(relative_c) + pi;
  %   phi3 = asin( relative_c ) - 2*pi*f
  %   % 20% = 0.2 corresponds to a phase shift of   -0.1pi
  %   %      -0.2   ~   +0.1pi
  % end
  % sinus_phi = sin(phi)
  % phi = 1;
  result = {t+1, 'prev_phase_shift', phi}; % + 0.25*pi
end


function result = chest_c4( model, trace, parameters, t )
  f     = trace(t+1).breathing_f.arg{1};
  range = trace(t+1).used_chest_range.arg{1};
  phi   = trace(t+1).prev_phase_shift.arg{1};

  dt    = model.parameters.default.dt;
  min   = model.parameters.default.min_chest_c;
  max   = model.parameters.default.max_chest_c;
  avg_chest_c = min + ((max - min) / 2);
  % f = dt * f;
  A = range / 2;
  % phi = phi * -1;
  curr_chest_c = A * sin(2*pi* f * dt + phi) + avg_chest_c;

  result = {t+1, 'chest_c4', curr_chest_c};
end





% starting direction can now be calculated using the new wave formula
% with the phase shift, but with t=2

function result = starting_dir( model, trace, parameters, t )
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
  if x1 < x2
    curr_dir2 = '1 in';
  else
    curr_dir2 = '3 out';
  end

  % old method
  %the 'observed' starting_dir    -> meetfout bij pieken
  %the starting directing of the breathing cycle  ()= current dir on t = t)
    %(this direction can change during dt)

  curr_chest_c = trace(t+1).chest_c4.arg{1};
  prev_chest_c = trace(t).chest_c4.arg{1};
  % prev_chest_c = trace(t).chest_c.arg{1};
  % if t>2
  %   curr_chest_c = trace(t+1).chest_c.arg{1};
  % else
  %   curr_chest_c = prev_chest_c;
  % end

  if prev_chest_c < curr_chest_c
    curr_dir = '1 in';
  else
    curr_dir = '3 out';
  end
  result = {t+1, 'starting_dir', {curr_dir2}};
end











% --------------------------------------------------------
%
%   ANALYSIS
%
% --------------------------------------------------------



function result = obs_chest_c( model, trace, parameters, t )

  chest_c = l2.getall(trace, t, 'chest_c', {NaN}).arg{1};
  % for c = l2.getall(trace, t, 'chest_c', {NaN})
  %   chest_c = c.arg{1};
   result = {t+1, 'observe', predicate('chest_c',chest_c)};
  % end
end

function result = bel_chest_c( model, trace, parameters, t )

  for c = l2.getall(trace, t+1, 'observe', predicate('chest_c', NaN))
    chest_c = c.arg{1}.arg{1};

    result = {t+1, 'belief', predicate('chest_c',chest_c)};
  end
end

function result = bel_chest_pos( model, trace, parameters, t )
  % for c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN))
  %   curr_chest_c = c.arg{1}.arg{1};
  %   for c2 = l2.getall(trace, t, 'belief', predicate('chest_c', NaN))
  %     prev_chest_c = c2.arg{1}.arg{1};
  curr_chest_c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
  prev_chest_c = l2.getall(trace, t, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};

  margin = model.parameters.default.margin;

  if curr_chest_c > prev_chest_c + margin
    chest_pos = '1 in';
  elseif curr_chest_c < prev_chest_c - margin
  	chest_pos = '3 out';
  else
  	chest_pos = '2 rest';
  end

  result = {t+1, 'belief', predicate('chest_pos',{chest_pos})};
  %   end
  % end
end

function result = bel_chest_trans( model, trace, parameters, t )
  %check whether a breathing cycle has started

  %oude manier (geeft betere error msg)
  % for c = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN))
  %   curr_pos = c.arg{1}.arg{1};
  %   for c2 = l2.getall(trace, t, 'belief', predicate('chest_pos', NaN))
  %     prev_pos = c2.arg{1}.arg{1};

  % chest_pos = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  prev_pos = l2.getall(trace, t, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  curr_pos = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
  %param om in- of uitademen te meten
  mode = model.parameters.default.chest_trans_mode;

  %check if there is a change,
  % when there is, check if the breathing cycle has restarted
  if strcmp(curr_pos,prev_pos) %no change
    chest_trans = false;
  elseif strcmp(curr_pos,'1 in') && mode == 1
    chest_trans = true;
  elseif strcmp(curr_pos,'3 out') && mode == 3
    chest_trans = true;
  else
    chest_trans = false;
  end

  result = {t+1, 'belief', predicate('chest_trans',chest_trans)};
end

function result = graph_bel_chest_pos( model, trace, parameters, t )

  chest_pos = l2.getall(trace, t+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};

  if strcmp(chest_pos,'1 in')
    chest_pos_number = 1;
  elseif strcmp(chest_pos,'2 rest')
    chest_pos_number = 2;
  else
    chest_pos_number = 3;
  end
  result = {t+1, 'graph_bel_chest_pos', chest_pos_number};
end

function result = graph_bel_chest_trans( model, trace, parameters, t )
  c = l2.getall(trace, t+1, 'belief', predicate('chest_trans', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_chest_trans', c};
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
      val = l2.getall(trace, a+1, 'belief', predicate('chest_trans', NaN)).arg{1}.arg{1};
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

function result = bel_breathing_acc( model, trace, parameters, t )
  % acceleration of the breathing f
  % breathing f starts with 0, therefore acc = 0
  breathing_f = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
  interval = model.parameters.default.acc_interval;
  margin = model.parameters.default.acc_margin;

  if interval > t-1, interval = t-1; end;
  end_t = t - interval;
  v = [];

  for a=t:-1:end_t
    val = l2.getall(trace, a, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
    count = t-a;              %teller: aantal loops
    v(count+1) = val;
    a = a-1;
  end
  %note that v is reversed, all the new (older t) values were appended to the vector

  avg = sum(v) / length(v);
  if length(v) == 0 disp('acc - div by 0'); end;
  if breathing_f > avg + margin
    breathing_acc = 1;
  elseif breathing_f < avg - margin
    breathing_acc = -1;
  else
    breathing_acc = 0;
  end

  result = {t+1, 'belief', predicate('breathing_acc',breathing_acc)};
end

function result = graph_bel_breathing_acc( model, trace, parameters, t )
  acc = l2.getall(trace, t+1, 'belief', predicate('breathing_acc', NaN)).arg{1}.arg{1};
  result = {t+1, 'graph_bel_breathing_acc', acc};
end

function result = bel_breathing_pattern( model, trace, parameters, t )

  interval = model.parameters.default.pattern_interval;

  if interval > t-1, interval = t-1; end;
  end_t = t - interval;
  v = [];

  for a=t:-1:end_t
    val = l2.getall(trace, a, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
    count = t-a;              %teller: aantal loops
    v(count+1) = val;
    a = a-1;
  end
  %note that v is reversed, all the new (older t) values were appended to the vector

  increases = v(v==1);
  decreases = v(v==-1);
  if ~isempty(increases) && ~isempty(increases)
    pattern = 'irregular';
  else
    pattern = 'regular';
  end

  result = {t+1, 'belief', predicate('breathing_pattern',{pattern})};
end

function result = graph_bel_breathing_pattern( model, trace, parameters, t )
  p = l2.getall(trace, t+1, 'belief', predicate('breathing_pattern', {NaN})).arg{1}.arg{1};
  result = {t+1, 'graph_bel_breathing_pattern', {p}};
end


% function result = adr6(trace, params, t)
%     result = {};
%
%     for loneliness_belief = l2.getall(trace, t, 'belief', predicate('feeling_of_loneliness', NaN))
%         belief = loneliness_belief.arg{1}.arg{1};

% for c = l2.getall(trace, t, 'observe', {'chest_c', NaN})

% for c = l2.getall(trace, t, 'chest_c', {NaN})
% for c = l2.getall(trace, t, 'observe', predicate('chest_c', NaN))
  % chest_c = c.arg{1}.arg{1};
  % for loneliness_belief = l2.getall(trace, t, 'belief', predicate('feeling_of_loneliness', NaN))
  %     belief = loneliness_belief.arg{1}.arg{1};

% function result = adr3(trace, params, t)
%     result = {};
%
%     for observation = l2.getall(trace, t, 'observation', {predicate('performance', NaN)})
%         belief = observation.arg{1};
%
%         result = {result{:} {t+1, 'belief', belief}};
%     end
% end

%oude

%function result = chest_pos( model, trace, parameters, t )
%
%    for c = l2.getall(trace, t+1, 'chest_pos', {NaN})
%        chest_pos1 = c.arg{1};
%
%        for c2 = l2.getall(trace, t, 'chest_pos', {NaN})
%            chest_pos2 = c2.arg{1};
%
%            if strcmp(chest_pos1,'in')
%                chest_pos_new = 'rest';
%            elseif strcmp(chest_pos1,'out')
%                chest_pos_new = 'rest';
%            elseif strcmp(chest_pos1,'rest')
%                if strcmp(chest_pos2,'in')
%                    chest_pos_new = 'out';
%                else
%                    chest_pos_new = 'in';
%                end
%            end
%
%        result = {t+2, 'chest_pos', {chest_pos_new}};
%        end
%    end
%end




% This rule is from the old model regarding the presence of food and the
% behaviour of an animal to this presence (Can be disregarded for the
% modded version).
%function result = ddr1( model, trace, parameters, t )
%    move = trace(t).observes_food_at_p2 & trace(t).observes_no_screen;
%    result = {t+1, predicate('goes_to_p2', move)};
%end

% This function reads input from the lightsensor and in turn switches the
% arduino LED on.
%function result = ddr2 (model, trace, parameters, t)
%    result = {};
%
%    % Obtain lightsensor value from the trace
%    for sensorPred = trace(t).lightvalue
%        sensor_value = sensorPred.arg{1};
%
%        % Initialize light_pred value
%        light_pred = '';
%        % Define the pin that contains the LED
%        pin = 'D4';
%
%        % If the value of the light sensor drops below a certain threshold
%        if sensor_value <= 1
%            display(sensor_value)               % Display the value actually obtained by the sensor (just to check)
%            light_pred = 'On';                  % Set the predicate for the LED to on
%            writeDigitalPin(model.controller, pin, 0);% Set the LED to the appropriate value (on)
%        else
%            display(sensor_value)
%            light_pred = 'Off';
%            writeDigitalPin(model.controller, pin, 1);
%        end
%
%        % Return the predicate that indicates the LED status to the trace
%        % (in order to keep track of what happened during the simulation)
%        result = {t+1, 'light', {light_pred}};
%    end
%end
