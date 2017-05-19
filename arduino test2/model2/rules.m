function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

function result = ddr12( model, trace, parameters, t )
    global n;
    if t == n-1
        bb = [];
        hrr = [];
        for i=1:t
            b2 = trace(i).b2.arg{1};
            hr = trace(i).hr2.arg{1};
            bb(i) = b2;
            hrr(i) = hr;
        end   
           
        disp('saving traces')
        ra = int2str(rand(1,1)*100);
        name1 = strcat('data/bb_v',ra,'.csv');
        name2 = strcat('data/hr_v',ra,'.csv');
        csvwrite(name1,bb);
        csvwrite(name2,hrr);
%        xlswrite('data/hee.xls',bb);
    end
    %csvwrite('pinda.csv', oi)
    
    result = {t+1, 'empty', true};
end

%ADD RULES BELOW

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
%     result = {};
%     % toc
%     %h = trace(t).hr2.arg{1};
%     %b = trace(t).b2.arg{1};
%     %h
%     %b
    
%     % Obtain lightsensor value from the trace
%     for sensorPred = trace(t).lightvalue
%         sensor_value = sensorPred.arg{1};
        
%         % Initialize light_pred value
%         light_pred = '';
%         % Define the pin that contains the LED
        
%         % If the value of the light sensor drops below a certain threshold
%         if sensor_value <= 1
%             %display(sensor_value)               % Display the value actually obtained by the sensor (just to check) 
%             light_pred = 'On';                  % Set the predicate for the LED to on
%             %writeDigitalPin(model.controller, pin, 0);% Set the LED to the appropriate value (on)
%         else
%             %display(sensor_value)
%             light_pred = 'Off';
%             %writeDigitalPin(model.controller, pin, 1);
%         end
        
%         % Return the predicate that indicates the LED status to the trace
%         % (in order to keep track of what happened during the simulation)
%         result = {t+1, 'light', {light_pred}};
%     end
% end