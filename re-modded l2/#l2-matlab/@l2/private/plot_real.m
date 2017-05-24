function [figs, heights] = plot_real( model, predicate )
%PLOT_REAL Plot a predicate with a single real value

    % 1 single figure of height 25mm
    figs = figure;
    heights = 25;
    set(figs,'Visible','off');
          
    % Get all the predicates
    yPreds = {model.trace.(predicate)};
    
    % Create the x vector
    x = [0:length(yPreds)];
    xLength = length(yPreds);
    addX = arrayfun(@(x) [x NaN x],x(2:end),'uniformOutput',false);
    x = [x(1) cell2mat(addX)];
    
    % Create y vectors (multiple if multiple values for same t exist)
    nrLines = max(cellfun(@(x) length(x), yPreds));    
    yValues = nan(nrLines, xLength);
    for i=1:xLength
        for j=1:length(yPreds{i})
            yPred = yPreds{i}(j);
            yValues(j,i) = yPred.arg{1};
        end
    end
    
    
    % Plot the figures
    hold on
    for i=1:nrLines
        y = yValues(i,:);
         % Transform values, such that it becomes a flat line for each timestep
        y = [y y(end)];
        addY = arrayfun(@(x) [x x x],y(1:end-1),'uniformOutput',false);
        y = [cell2mat(addY) y(end)];

        % Plot the graph and add title as ylabel
        plot(x,y);         
    end
    hold off
    ylabel(predicate);
end

