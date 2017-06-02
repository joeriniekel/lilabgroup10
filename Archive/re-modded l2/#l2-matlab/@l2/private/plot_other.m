function [figs, heights] = plot_other( model, predicate )
%PLOT_OTHER Plot a predicate with a single (string) value
% For each unique value in predicate, a boolean like graph is returned

    % Create figure and size arrays 
    figs = [];
    heights = [];

    % Get all the unique values in the predicate
    values = cellfun(@(x) pred2str(x), {model.trace.(predicate)}, 'UniformOutput', false);
    uniqueValues = unique([values{:}]);
    if ~iscell(uniqueValues)
        uniqueValues = unique({values{:}});
    end
     
    % Go through each value and create plot
    for i=1:length(uniqueValues)
        
        % 1 single figure of 5mm
        figs = [figs figure];
        heights = [heights 5];
        set(figs(end),'Visible','off');
        
        % Find all the values for this unique value
        y = cellfun(@(x) any(strcmp(x, uniqueValues{i})), values);

        % Create the graph
        hold on
        [xxTrue,yyTrue] = stairs([0:length(y)],[y 0]);
        patch([0 xxTrue']',[0 yyTrue']',[0 0 1], 'LineStyle','none');
        hold off


        % Fix the y-axis and add label
        set(gca,'YLim', [-0.5 1.5]);
        set(gca,'YTick', [-0.5 1.5]);
        set(gca,'YTickLabel',[]);
        set(gca, 'Layer','top');
        ylabel(uniqueValues{i})
    end
end

