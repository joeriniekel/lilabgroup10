function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

function result = anxiety( model, trace, parameters, t )

    for sf = l2.getall(trace, t, 'sitfac', {NaN})
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

    %for p = l2.getall(trace, t+1, 'physical_state', {NaN})
    %    ps = p.arg{1};

    anxiety = trace(t+1).anxiety.arg{1};
    bhr = model.parameters.default.bhr;
    ps = trace(t).physical_state.arg{1}; % t+1 doesn't work with this syntax + scenario values
    a = model.parameters.default.anxiety_hr;

    hr_new = (bhr * ps) + (a * anxiety);
    %disp('hr')
    %disp(hr_new)
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

  for c = l2.getall(trace, t+1, 'chest_c', {NaN})
    current_chest_c = c.arg{1};
      for c2 = l2.getall(trace, t, 'chest_c', {NaN})
        prev_chest_c = c2.arg{1};

        margin = model.parameters.default.margin;

      	if current_chest_c > prev_chest_c + margin
          chest_pos = '1 in';
        elseif current_chest_c < prev_chest_c - margin
        	chest_pos = '3 out';
        else
        	chest_pos = '2 rest';
        end

        result = {t+2, 'chest_pos', {chest_pos}};
      end
    end
end




% --------------------------------------------------------
%
%   ANALYSIS
%
% --------------------------------------------------------

function result = obs_chest_c( model, trace, parameters, t )

  for c = l2.getall(trace, t+1, 'chest_c', {NaN})
    chest_c = c.arg{1};

    result = {result{:} {t+1, 'observe', predicate('chest_c', chest_c)}};
  end
end



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
