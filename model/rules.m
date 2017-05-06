function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

function result = anxiety( model, trace, parameters, t )

  sitfac2 = trace(t).sitfac.arg{1}; %gives errors
        

    for sf = l2.getall(trace, t, 'sitfac', {NaN}) %t+1? todo
        sitfac = sf.arg{1};

        prev_anxiety = trace(t).anxiety.arg{1};
        %sitfac = trace(t).sitfac.arg{1}; %gives errors
        breathing_f = trace(t).breathing_f.arg{1};

        decay = model.parameters.default.anxiety_decay;
        disfac = model.parameters.default.disfac;
        s = model.parameters.default.sitfac_anxiety;
        b = model.parameters.default.breathing_anxiety;

        %dispositional factors affect the intensity of anxiety provoking
        %situations
        new_anxiety = prev_anxiety * decay + ...
            disfac * (s * sitfac + b * breathing_f);

        result = {t+1, 'anxiety', new_anxiety};
    end
end

function result = hr( model, trace, parameters, t )

  anxiety = trace(t+1).anxiety.arg{1};
  bhr = model.parameters.default.bhr;
  ps = trace(t).physical_state.arg{1}; % t+1 doesn't work with this syntax + scenario values
  a = model.parameters.default.anxiety_hr;

  hr_new = (bhr * ps) + (a * anxiety);

  result = {t+1, 'hr', hr_new};
end

function result = breathing_f( model, trace, parameters, t )

    %breathing_f = trace(t+1).breathing_f.arg{1};
    hr = trace(t+1).hr.arg{1};
    h = model.parameters.default.hr_breathing;

    breathing_f_in_bpm = h * hr;
    breathing_f_new = breathing_f_in_bpm /60;
    %disp('breathing f in bpm')
    %disp(breathing_f_in_bpm)
    result = {t+1, 'breathing_f', breathing_f_new};
end

function result = chest_c( model, trace, parameters, t )

    dt = model.parameters.default.dt;
    max_c = model.parameters.default.max_chest_c;
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

    result = {t+1, 'chest_c', chest_c};
end

function result = chest_pos( model, trace, parameters, t )
  % anxiety = trace(t).chest_c.arg{1};
  % for c = l2.getall(trace, t+1, 'chest_c', {NaN})
  %   current_chest_c = c.arg{1};
  %     for c2 = l2.getall(trace, t, 'chest_c', {NaN})
  %       prev_chest_c = c2.arg{1};
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
    %   end
    % end
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

% function result = chest_trans_array( model, trace, parameters, t )
%
%   % make a vector with old values
%   % start = 1;
%   % max = 10;
%   % if t > 10
%   %   start = t-max;
%   % end
%   % vector = [];
%   % cell = {};
%   % for i=start:t
%   %   val = l2.getall(trace, i+1, 'belief', predicate('chest_c', NaN)).arg{1}.arg{1};
%   %   val2 = l2.getall(trace, i+1, 'belief', predicate('chest_pos', NaN)).arg{1}.arg{1};
%   %   cell{i} = [val2];
%   % end
%   chest_trans_array = trace(t-x).chest_trans_array.arg{1};
%   % c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN))
%   % disp(c)
%
%   % for c = l2.getall(trace, t+1, 'belief', predicate('chest_c', NaN))
%   % disp(chest_c_array);
%
%   result = {t+1, 'chest_trans_array', chest_trans_array};
%   %   result = {t+1, 'chest_c_array', chest_c_array};
%   % end
% end

function result = bel_breathing_f( model, trace, parameters, t )

  % chest_pos = l2.getall(trace, t+1, 'belief', predicate('chest_trans', NaN)).arg{1}.arg{1};
  % start = 1;
  % max = 10;
  % if t > 10
  %   start = t-max;
  % end
  % v = [];
  % % cell = {};
  % for i=start:t
  %   val = l2.getall(trace, i+1, 'belief', predicate('chest_trans', NaN)).arg{1}.arg{1};
  %   % cell{i} = [val];
  %   v(i) = val;
  % end
  %raw calculation (incl. incomplete breathing cycles)
  % n_cycles = sum(v);  % count transitions to 'in'; number of instances that are true in the vector
  % breathing_f = n_cycles / max;
  %
  % first_start_of_cycle = find(v==1,1); %returns the index number

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
    %note that v is reversed, all the new (older t) values were appenden to the vector

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
  interval = 5;
  margin = 10;
  if interval > t -1, interval = t -1; end;
  end_t = t - interval;
  v = [];
  for a=t:-1:end_t
    val = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN)).arg{1}.arg{1};
    count = t-a;              %teller: aantal loops
    v(count+1) = val;
    a = a-1;
  end
  %note that v is reversed, all the new (older t) values were appenden to the vector
  % v
  % disp('sum')
  % sum(v)
  % length(v)
  avg = sum(v) / length(v);
  if length(v) == 0 disp('acc - div by 0'); end;
  breathing_f
  if breathing_f > avg + margin
    breathing_acc = 'increasing';
  elseif breathing_f < avg - margin
    breathing_acc = 'decreasing';
  else
    breathing_acc = 'stable';
  end

  breathing_acc = 5;
  result = {t+1, 'belief', predicate('breathing_acc', breathing_acc)};
end



% function result = bel_breathing_acc( model, trace, parameters, t )

%   chest_c = l2.getall(trace, t, 'breathing_f', NaN).arg{1};
%   for x  = l2.getall(trace, t+1, 'belief', predicate('breathing_f', NaN));
%     breathing_f = x.arg{1}.arg{1};
%     breathing_f


%     breathing_acc = 3;
%   % for c = l2.getall(trace, t, 'chest_c', {NaN})
%   %   chest_c = c.arg{1};
%     result = {t+1, 'belief', predicate('breathing_acc', breathing_acc)};

%     % result = {t+1, 'observe', predicate('chest_c',chest_c)};
%   end
% end


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

%oude cyclus
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
