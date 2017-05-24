function [figs, heights] = plot_boolean( model, pred )
%PLOT_BOOLEAN Plot a graph for a single boolean value

    % 1 single figure of 5mm
    figs = figure;
    heights = 5;
    set(figs,'Visible','off');
    
    % Get all the values and make false -.5 and true 0.5    
    yPreds = {model.trace.(pred)};
    xLength = length(yPreds);
     
    % yValues, true at yValues(1,:), false at yValues(2,:)
    yFalse = zeros(1, xLength);
    yTrue = zeros(1, xLength);
    for i=1:xLength
        for j=1:length(yPreds{i})
            yPred = yPreds{i}(j);
            if yPred.arg{1}
                yTrue(i) = 1;
            elseif ~yPred.arg{1}
                yFalse(i) = -1;
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
    ylabel(pred);
    
end