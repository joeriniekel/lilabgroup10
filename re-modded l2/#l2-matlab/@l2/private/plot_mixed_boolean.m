function [figs, heights] = plot_mixed_boolean( model, pred )
%PLOT_MIXED_BOOLEAN Plot a boolean graph for each unique combination of
%preceding arguments, boolean value is always last

    
    % Get the sets to plot, that is those unique values ignoring the BOOLEAN
    pAll = pred2str([model.trace.(pred)]);
    sets = unique(regexprep(pAll, ', (false|true)}$', '}'), 'rows');
    if ~iscell(sets)
        sets = {sets};
    end
    
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
            setIndex = find(strcmp(regexprep(pred2str(yPred), ', (false|true)}$', '}'), sets));
            ySorted{setIndex}{i} = [ySorted{setIndex}{i} yPred];
        end
    end

    % Create figure and size arrays 
    figs = [];
    heights = [];
    
    % Go through each set and create plot
    for i=1:length(sets)
        
        % 1 single figure of 5mm
        figs = [figs figure];
        heights = [heights 5];
        set(figs(end),'Visible','off');
        
        % Get all the values and make false -.5 and true 0.5    
        yPreds = ySorted{i};
        xLength = length(yPreds);

        % yValues, true at yValues(1,:), false at yValues(2,:)
        yFalse = zeros(1, xLength);
        yTrue = zeros(1, xLength);
        for j=1:xLength
            for k=1:length(yPreds{j})
                yPred = yPreds{j}(k);
                if yPred.arg{end}
                    yTrue(j) = 1;
                elseif ~yPred.arg{end}
                    yFalse(j) = -1;
                end
            end
        end

        % Plot the figures
        hold on
        [xxFalse,yyFalse] = stairs([0:length(yFalse)],[yFalse 0]);
        [xxTrue,yyTrue] = stairs([0:length(yTrue)],[yTrue 0]);
        patch([0 xxFalse']',[0 yyFalse']',[0.7 0.85 0.9]);
        patch([0 xxTrue']',[0 yyTrue']',[0 0 1]);
        hold off

        % Fix the y-axis and add label
        set(gca,'YLim', [-1.5 1.5]);
        set(gca,'YTick', [-1 1]);
        %set(gca,'YTickLabel',{'false', 'true'});
        set(gca,'YTickLabel',[]);
        set(gca, 'Layer','top');
        ylabel(sets{i});
    end

end