function [ score ] = fitness( model )
%ESTIMATOR Output a parameter score based on the trace

    % Read the results to compare trace to
    testData = csvread([mfilename('fullpath') '.csv'],1,0);
    
    % Make sure that both the trace and desired results have an equal length
    if length([model.trace.erl]) >= length(testData(:,1)) 
        target = testData(:,1);
        simulated = [model.trace(1:length(testData(:,1))).erl]';
    else
        target = testData(1:length([model.trace.erl]),1);
        simulated = [model.trace.erl]';
    end
    
    % Fitness is a simple mean squared error
    score = sqrt( mean((target-simulated).^2) );
    
    % If for some reason no score can be calculated, fitness is as high as
    % possible (i.e. the worst result possible)
    if isnan(score) 
        score = intmax;
    end
end

