function [ps, height] = plot_panel( model, toPlot, target )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Set textinterpreter for plotting
    set(0, 'DefaulttextInterpreter', 'none');
    
    % Check if any of the to plot predicates have no value
    emptyPreds = false(1,length(toPlot));
    for i=1:length(toPlot)
        emptyPreds(i) = isempty([model.trace.(toPlot{i})]);
    end
    toPlot(emptyPreds) = [];
    nrPlots = length(toPlot);
    
    %% Get plotting functions
    funcPlot = cell(nrPlots,1);
    % Temporarily change dir to model folder
    prevDir = pwd;
    clean__ = onCleanup(@(~)cd(prevDir));
    cd(model.model);
    for i=1:nrPlots
        
        % If a function of name plot_<predicate> exists, use that
        funcName = ['plot_' toPlot{i}];
        if exist(funcName)
            funcPlot{i} = str2func(funcName);
            
        % Otherwise determine function based on sort
        else
            predSort = model.getSort(toPlot{i});
            % cell array
            if iscell(predSort)
                funcPlot{i} = @plot_cell;
            % REAL
            elseif strcmp(predSort, 'REAL') 
                funcPlot{i} = @plot_real;
            % BOOLEAN
            elseif strcmp(predSort, 'BOOLEAN')
                funcPlot{i} = @plot_boolean; 
            % OTHER (string predicate)
            else
                funcPlot{i} = @plot_other; 
            end
        end
    end
    delete(clean__);
    
    %% Plot the graphs
    
    % Setup some variables
    xLim = length(model.trace);
    rowMargin = 2;
    pMargin = [2*rowMargin 2*rowMargin 2*rowMargin 3*rowMargin];
    height = pMargin(2) + pMargin(4);
    onScreen = strcmp(target, 'Screen');
    
    % Create a panel
    p = newPanel();
    ps = {};
    p.margin = pMargin;
    newWindow = false;
    
    wb=waitbar(0,'Plotting..');
    clean__ = onCleanup(@()close(wb));
    % Add all plots one by one
    for i=1:length(funcPlot)
        % Get figures and sizes from plotting function
        [figs, heights] = funcPlot{i}(model, toPlot{i});

        % Copy figures to the panel and set labels, ticks, etc.
        for j=1:length(figs)
            p.pack({{heights(j)}});
            pI = length(p.children);
            p(pI).margin = rowMargin;
            xAxis = '';
            
            % Get the position
            pPos = p(pI).position;

            %Save the top margin if its the first plot
            if pI==1
                margintop = 1-pPos(2)-pPos(4);
            end

            % Check if it still fits
            if onScreen & pPos(2) <= margintop
                newWindow = true;
            elseif ~onScreen & (height(end) + heights(j) + rowMargin + margintop) > 250
                newWindow = true;
            end
            if newWindow
                % Delete the last graph
                delete(p(pI));

                % Add bottom x-axis to previous graph
                addAxes(p, pI-1, 'bottom');

                % Create new panel and add current graph
                ps = [ps {p}];
                p = newPanel();
                p.margin = pMargin;
                p.pack({{heights(j)}});
                pI = length(p.children);
                p(pI).margin = rowMargin;
                xAxis = 'top';
                height = [height pMargin(2) + pMargin(4)];
                newWindow = false;
            end
            
            % Save the height and add axis if required
            height(end) = height(end) + heights(j) + rowMargin;
            if i==length(funcPlot) && j==length(figs)
                xAxis = 'bottom';
            elseif i==1 && j==1
                xAxis = 'top';
            end
            
            % Copy the graph
            ml = copyFig(figs(j), p, pI, xLim);
            addAxes(p, pI, xAxis);
            pMargin(1) = max(pMargin(1), ceil(abs(ml*10)));
            
            % Close the old figure
            close(figs(j));
        end
        waitbar(i/length(funcPlot));
    end
    % Add panel to array
    ps = [ps {p}];
    
    delete(clean__);
    
    % Finalize panel(s)
    for i=1:length(ps)
        ps{i}.margin = pMargin;
    end
    set(findobj(ps{i}.figure,'Type','axes','Tag','legend'), 'Orientation','Horizontal', 'Location','best', 'FontSize', 8);

end

function marginleft = copyFig(fig, p, i, xLim)
    
    % Parse fig
    aFig = [];
    aLeg = [];
    fc = get(fig,'Children');
    for j = length(fc):-1:1
        if strcmp(get(fc(j),'Type'),'axes') && ~strcmp(get(fc(j),'Tag'),'legend')
            aFig = [aFig fc(j)];
        elseif strcmp(get(fc(j),'Type'),'axes') && strcmp(get(fc(j),'Tag'),'legend')
            aLeg = fc(j);
        end
    end

    % Copy figure and set axis
    p(i).select(aFig);
    set(gca,'Box', 'on');        
    set(gca,'XLim', [0 xLim]);
    set(gca,'XTickMode','auto',...
        'XTickLabel',[],...
        'xgrid','on',...
        'FontSize', 8);
    set(gca,'ticklength',0*get(gca,'ticklength'));
    set(get(gca,'YLabel'), ...
        'Rotation', 0,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle', ...
        'FontSize', 8, 'FontWeight', 'bold')
    set(get(gca, 'YLabel'), 'Unit', 'centimeters');
    labelPos = get(get(gca, 'YLabel'), 'Extent');
    marginleft = labelPos(1);
    
    % Copy colormap
    colormap(colormap(fig));
        
    % Check for legend and copy if exists 
    if ~isempty(aLeg)
        lData = get(aLeg, 'Userdata');
        legend(lData.handles, lData.lstrings);
    end
end

function addAxes(p,i,loc)
    p(i).select();
    if strcmp(loc, 'top')
        set (gca, 'XAxisLocation', 'top', 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'FontSize', 8)
    elseif strcmp(loc, 'bottom')
        set(gca, 'XAxisLocation', 'bottom', 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'FontSize', 8)
    elseif strcmp(loc, 'both')
        set(gca, 'XAxisLocation', 'bottom', 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'FontSize', 8)
        % TODO: copy x axis to top as well?
    end
end

function p = newPanel() 
    
    % Create the figure
    figure;
    p = panel('no-manage-font');
    
    % Set position
    movegui(p.figure(), 'northwest');
    scnSize  = get(0,'ScreenSize');
    set(p.figure(),'outerposition',[scnSize(1) scnSize(2)+50 min(900, scnSize(3)) scnSize(4)-50]);
end