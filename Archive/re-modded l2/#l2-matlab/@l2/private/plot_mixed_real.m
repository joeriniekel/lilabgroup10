function [figs, heights] = plot_mixed_real( model, pred )
%PLOT_MIXED_REAL Plot a predicate with two values one of which is a REAL

    % Get the index of REAL
    sorts = model.predicates.(pred);
    iReal = 0;
    iOther = 0;
    for i=1:length(sorts)
        if strcmp(sorts{i}, 'REAL')
            iReal = i;
        else
            iOther = i;
        end
    end
    % Get the sets to plot, that is those unique values ignoring the REAL
    pAll = [model.trace.(pred)];
    sets = unique(cellfun(@(x) x{1}, {pAll.arg}, 'UniformOutput', false));
   
    % Get the values and divide the for the different sets
    yAll = {model.trace.(pred)};
    ySorted = cell(1, length(sets));
    
    % Initialise each set
    for i=1:length(sets)
        ySorted{i} = cell(1,length(yAll));
    end
    
    % Go through each timepoint
    for i=1:length(yAll)
        % Go through each predicate
        for j=1:length(yAll{i})
            % Add the value to the corresponding set
            yPred = yAll{i}(j);
            setIndex = find(strcmp(yPred.arg{iOther}, sets));
            ySorted{setIndex}{i} = [ySorted{setIndex}{i} yPred];
        end
    end

    % 1 single figure, height depends on number of sets
    figs = figure;
    heights = 20 + 5*length(sorts);
    set(figs,'Visible','off');
    ColOrd = get(gca,'ColorOrder');
    [m,n] = size(ColOrd);
    legendLines = [];
    
    % Create the x vector
    x = [0:length(yAll)];
    xLength = length(yAll);
    addX = arrayfun(@(x) [x NaN x],x(2:end),'uniformOutput',false);
    x = [x(1) cell2mat(addX)];
    
    % Now plot each set as it was a single real value in a single figure
    hold all
    for s=1:length(sets)
        yPreds = ySorted{s};
        ColRow = rem(s,m);
        if ColRow == 0
          ColRow = m;
        end
        
        % Create y vectors (multiple if multiple values for same t exist)
        nrLines = max(cellfun(@(x) length(x), yPreds));    
        yValues = nan(nrLines, xLength);
        for i=1:xLength
            for j=1:length(yPreds{i})
                yPred = yPreds{i}(j);
                yValues(j,i) = yPred.arg{iReal};
            end
        end
        
        % Plot the figures
        for i=1:nrLines
            y = yValues(i,:);
             % Transform values, such that it becomes a flat line for each timestep
            y = [y y(end)];
            addY = arrayfun(@(x) [x x x],y(1:end-1),'uniformOutput',false);
            y = [cell2mat(addY) y(end)];

            % Plot the graph and add title as ylabel
            if i==1
                legendLines = [legendLines plot(x,y, 'Color', ColOrd(ColRow,:))];
            else
                plot(x,y, 'Color', ColOrd(ColRow,:))
            end
        end
        
    end
    
    % Add a label and legend
    ylabel([pred '{' sorts{iOther} '}']);
    legend(legendLines, sets{:}, 'Orientation','Horizontal', 'Location','Best')
    hold off
end

