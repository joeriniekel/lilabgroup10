function trace = step( model, trace, params, rules, t, scnName)
%STEP Perform one iteration of model
%   Detailed explanation goes here
    
    %get function handles from rules file
    fcns = feval(rules);
    
    %update trace with data from current timepoint t
    trace = setScenarioLoop( model, trace, scnName, t);
    
    %go through each rule
    for i=1:length(fcns)
        fHandle =fcns{i};
        try
            result = fHandle( model, trace, params, t );
        catch err
            causeErr = MException('L2:Step', 'Code from rule ''%s'' (line %i).', func2str(fHandle), err.stack(1).line);
            err = addCause(err, causeErr);
            throwAsCaller(err);
        end
        
        try
            trace = l2add(model, trace, result);
        catch err
            causeErr = MException('L2:Step', 'Result from rule ''%s''.', func2str(fHandle));
            err = addCause(err, causeErr);
            throwAsCaller(err);
        end
    end
end