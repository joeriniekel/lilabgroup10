function result = pred2str( pred )
%PRED2STR Convert predicate format to string
%   Convert a given predicate to string. Input is a single predicate, in
%   which case a string is returned, or an array of predicates in which
%   case a cell array of strings is returned.
%
%   See also PREDICATE, ISPREDICATE, PREDCMP, STR2PRED

    if length(pred) == 1
        result = predicate2string(pred);
    else
        result = cell(1, length(pred));
        for i=1:length(pred)
            result{i} = predicate2string(pred(i));
        end
    end
end

function result = predicate2string(pred)
    name = pred.name;
    values = values2string(pred.arg);
    result = [name '{' values '}'];
end

function result = values2string(arg)
    strings = cell(1, length(arg));
    for i=1:length(arg)
        if isa(arg{i}, 'predicate')
            strings{i} = predicate2string(arg{i});
        else
            strings{i} = value2string(arg{i});
        end
    end
    result = strjoin(strings, ', ');
end

function result = value2string( value )
%value2string Convert single value to a string
    %single value (logical)
   	if islogical(value)
        if value
            result = 'true';
        else
            result = 'false';
        end
        
    %single value (NaN)
    elseif isnan(value)
        result = '~';
    
    %single value (numerical)    
    elseif isnumeric(value)
        result = num2str(value);
        
    %single value (string)
    else
        result = value;
    end
end