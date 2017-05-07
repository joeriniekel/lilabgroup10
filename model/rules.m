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
  tijdeleijke_param = 0;

  anxiety_regulation = disfac + b * (1 - breathing_f + tijdeleijke_param);

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

  anxiety = trace(t+1).anxiety.arg{1};
  bhr = model.parameters.default.bhr;
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
  breathing_f = breathing_f_in_bpm / 60;

  result = {t+1, 'breathing_f', breathing_f};
end


% chest_c veroorzaakt nu schokken
% een betere modulatie zou zijn:
%   bepaal huidige transition:  in- of uitademen
%   bepaal hr, bepaal breathing f
%   bepaal de snelheid waarmee de chest pos moet bewegen
%     (f/afstand)
%   bepaal aan de hand daarvan de nieuwe chest pos (e.g. oude pos + 2cm)
%
%   op deze manier is het modeleren van afwijkingen ook logischer

% hr -> breathing f
% breathing f -> breathing intensity -> used chest c range
% breathing f + intensity -> chest movement speed (average)



%breathing_intensity
% the relation of this concept is unknown. it should look like this:

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

% function result = chest_speed( model, trace, parameters, t )
%   breathing_f = trace(t+1).breathing_f.arg{1};
%   used_chest_range = trace(t+1).used_chest_range.arg{1};
%   chest_speed = used_chest_range / breathing_f;
%   result = {t+1, 'chest_speed', chest_speed};
% end
% chest_speed?

function result = chest_pos_phi( model, trace, parameters, t )
  prev_chest_c = trace(t).chest_c.arg{1};
  prev_chest_pos = trace(t).chest_pos.arg{1};
  f = trace(t+1).breathing_f.arg{1};
  chest_range = trace(t+1).used_chest_range.arg{1};
  dt = model.parameters.default.dt;
  min = model.parameters.default.min_chest_c;
  max = model.parameters.default.max_chest_c;

  avg_chest_c = min + (max - min) / 2;    % the 'zero-line' of the wave/cycle (x(t) = 0)
  A = chest_range / 2;
  margin_above = avg_chest_c + A;
  margin_below = avg_chest_c - A;

  %the formula for a standard waveform is
  % x(t) = A * sin(2*pi*f*t + phi)
  %   phi in range[0,2*pi]
  % thus:
  % chest_c(t) = avg_chest_c + A * sin(2*pi*f*t + phi)
  % new_chest_c = avg_chest_c + A * sin(2*pi*f*dt + phi)

  % calculate phi of the prev chest_c diverted of the new breathing f
  % x(t) = A * sin(2*pi*f*t + phi)
  % x(0) = A * sin(2*pi*f*0 + phi)
  % x(0) = A * sin(0 + phi)
  % x(0) = A * sin(phi)
  % let x(0) = prev_x
  % prev_x = A * sin(phi)
  % prev_x / A = sin(phi)
  % phi = asin(prev_x / A)
  prev_x = (prev_chest_c - avg_chest_c) / chest_range
    % in range [-1,1]
  % a value of 1 corresponds to phi = 0.5pi or phi = 1.5pi
  % (this depends on whether the agent is breathing in our out)
  phi = asin( prev_x / A )
  % prev_chest_c
  new_chest_c = avg_chest_c + A * sin(2*pi*f*dt + phi)

  % phi in range[-0.5; 0.5]
  % phase = positive number
  % reduced/instantaneous phase = in range[0,1]
  % phase = mod(n,1);     %the phase of the cycle: a number between 0 and 1

  % phase of breathing cycle   chest_c
  % 0 and 0.5                  'zero-line' = avg_chest_c
  % 0.25                       max
  % 0.75                       min

  %calculate the phase of the prev chest_c diverted of the new breathing f
  % prev_chest_pos
  % prev_chest_c
  % avg_chest_c
  % amplitude
  % pi
  % phase = calc_phase(prev_chest_pos,prev_chest_c,avg_chest_c,amplitude)
  % if strcmp(prev_chest_pos, '1 in') && prev_chest_c > avg_chest_c
  %   %phase is between 0 and 0.25
  %   relative_chest_c = (prev_chest_c - avg_chest_c) / amplitude;
  %   phase = relative_chest_c / 4;
  % elseif strcmp(prev_chest_pos, '1 in') && prev_chest_c < avg_chest_c
  %   %phase is between 0.75 and 0
  %   relative_chest_c = (avg_chest_c - prev_chest_c) / amplitude;
  %   phase = (relative_chest_c / 4) + 0.75;
  % elseif strcmp(prev_chest_pos, '3 out') && prev_chest_c > avg_chest_c
  %   %phase is between 0.25 and 0.5
  %   relative_chest_c = (prev_chest_c - avg_chest_c) / amplitude;
  %   phase = (relative_chest_c / 4) + 0.25;
  % else
  %   % strcmp(prev_chest_pos, '3 out') && prev_chest_c < avg_chest_c
  %   %phase is between 0.25 and 0.5
  %   relative_chest_c = (avg_chest_c - prev_chest_c) / amplitude;
  %   phase = (relative_chest_c / 4) + 0.5;
  % end

  %the phase of the previous chest_c is now calculated
  %now calculate the new chest_c
  % new_phase = phase + chest_speed / dt;
  % new_chest_c = phase
  %if chest_c > margin_above
  % change direction
  %if chest_c < margin_below
  % change direction


  %option: use sinus waveform instead of triangle
  % chest_pos = '1 in'; % POS;          {1 in, 2 rest, 3 out}
  result = {t+1, 'chest_phi', phi};
end

function result = chest_c( model, trace, parameters, t )
  prev_chest_c = trace(t).chest_c.arg{1};
  f = trace(t+1).breathing_f.arg{1};
  chest_range = trace(t+1).used_chest_range.arg{1};
  dt = model.parameters.default.dt;
  min = model.parameters.default.min_chest_c;
  max = model.parameters.default.max_chest_c;

  avg_chest_c = min + (max - min) / 2;    % the 'zero-line' of the wave/cycle (x(t) = 0)
  A = chest_range / 2;

  % prev_x = (prev_chest_c - avg_chest_c) / chest_range;
  %
  % phi = asin( prev_x / A );
  phi = trace(t+1).chest_phi.arg{1};
  new_chest_c = avg_chest_c + A * sin(2*pi*f*dt + phi);
  % chest_c = 70;

  result = {t+1, 'chest_c', new_chest_c};
end

% oude chest_c
% function result = chest_c( model, trace, parameters, t )
%
%     dt = model.parameters.default.dt;
%     max_c = model.parameters.default.max_chest_c;
%     min_c = model.parameters.default.min_chest_c;
%     breathing_f = trace(t+1).breathing_f.arg{1};
%     %breathing_f = n/s
%     %n = breathing_f * s;
%
%     n = breathing_f * t * dt;   %number of breathing cycles
%     phase = mod(n,1);     %the phase of the cycle: a number between 0 and 1
%     %in this case the breathing cycle always starts at t=0
%
%     if phase < 0.5
%         phi = phase;
%     else
%         phi = 1 - phase;
%     end
%     chest_c = max_c * 2 * phi;
%
%     result = {t+1, 'chest_c', chest_c};
% end

% oude chest_pos
function result = chest_pos( model, trace, parameters, t )
  curr_chest_c = trace(t+1).chest_c.arg{1}
  prev_chest_c = trace(t).chest_c.arg{1}
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
      disp('no cycles found')
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
