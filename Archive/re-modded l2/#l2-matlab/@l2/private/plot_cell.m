function [figs, heights] = plot_cell( model, predicate )
%PLOT_CELL Determine how to plot a predicate with multiple values    

    sorts = model.getSort(predicate);

    
    if strcmp(sorts{end}, 'BOOLEAN')
    % plot boolean like graph for each unique combination    
        [figs, heights] = plot_mixed_boolean( model, predicate );
        
    elseif length(sorts) == 2 ...
            && sum(cellfun(@(x) strcmp(x, 'REAL'), sorts)) == 1 ...
            && length(unique(pred2str([model.trace.(predicate)]))) > 4
    % predicate with a single REAL value and multiple results
        [figs, heights] = plot_mixed_real( model, predicate );
        
    else
    % plot anything else as other
        [figs, heights] = plot_other( model, predicate );
    end
    
    
    % TODO, check for a nested single REAL

end

