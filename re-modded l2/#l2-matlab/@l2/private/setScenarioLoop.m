function trace = setScenarioLoop( model, trace, scnName , t)
%SETSCENARIOLOOP Load new scenario entries into the trace for every new
%point in time. This function should be called from the step function.

    scnToSet = model.scenarios.(scnName);
    
    % for all fields
    fields = fieldnames(scnToSet);
    for i=1:size(fields)
        % add all nonempty rows
        for j=1:length(scnToSet.(fields{i}))
            interval = scnToSet.(fields{i}){j}{1};
            pred = scnToSet.(fields{i}){j}{2};
            % if interval contains the current time point set it to this
            % time point
            currentTime = strfind(interval, t);
            
            %if predicate value contains 'rtinput_sensorfilename' retrieve
            %the predicate using sensor data from the file
            if strfind(pred, 'rtinput_')
                sensorPred = strsplit(pred.arg{1}, '_'); % split off the filename from the rtinput_sensorfilename
                sensorHandle = str2func(sensorPred{2}); % get the sensorfilename from the predicate and create a function handle out of it.
                
                cd(model.model) % change directory to the model folder to access the sensor script
                try
                    pred = sensorHandle(model.controller);
                    cd('..') % change directory back after successfully updating the predicate value using the sensor
                catch err
                    cd('..') % change the directory back in case of a failure
                    causeErr = MException('L2:SetScenario', 'Unable to invoke result of sensor script ''%s''.', sensorPred{2});
                    err = addCause(err, causeErr);
                    throwAsCaller(err)
                end
            end
            
            % if interval contains Inf, replace with the current time point
            indexInf = strfind(interval, 'Inf');

            %%%
            % if t == 57
            %     disp('t is 57')
            %     disp(pred)
            %     interval
            % end
            %%%

            if indexInf > 0 
                split_interval = strsplit(interval,':');
                t_start = strsplit(split_interval{1}, '[');
                t_start = str2double(t_start(2));
                if t >= t_start
                    interval = t;
                end
            end
            if currentTime > 0
                interval = t;
            end
            try
                trace = l2add(model, trace, {interval, pred});
            catch err
                causeErr = MException('L2:SetScenario', 'Entry in scenario ''%s''.', scnName);
                err = addCause(err, causeErr);
                throwAsCaller(err);
            end
        end
    end
    
end