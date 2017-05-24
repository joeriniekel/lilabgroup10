function result = predcmp( pred1, pred2 )
%PREDCMP Compare 2 predicates on equality
%   Compare whether 2 predicates are equal. Input can be either in string
%   or predicate form. In a string, locations with a '~' are ignored, while
%   in a predicate values to be ignored are represented by a 'NaN'.
%
%   See also PREDICATE, ISPREDICATE, PRED2STR, STR2PRED
    try
        if ischar(pred1)
            pred1 = str2pred(pred1);
        end
        if ischar(pred2)
            pred2 = str2pred(pred2);
        end
    catch err
        %throwAsCaller(err)
        rethrow(err)
    end

    result = comparepredicate(pred1, pred2);
end


function result = comparepredicate(pred1, pred2 )
    % compare name and values on equality
    name = comparename(pred1.name, pred2.name);
    values = comparevalues(pred1.arg, pred2.arg);
    % result is logical combination of both
    result = name && values;
end

function result = comparename(name1, name2)
    %either value is not relevant
    if eithernan( name1, name2 )
        result = true;
        return
    
    % compare the two name strings
    else
        result = strcmp(name1, name2);
    end
end

function result = comparevalues(arg1, arg2)
    % values need to have equal number to be the same
    if length(arg1) ~= length(arg2)
        result = false;
        return
    end

    % compare the values one by one until a mismatch is found
    for i=1:length(arg1)
        if eithernan( arg1{i}, arg2{i} )
            continue
        end
        
        if isa(arg1{i}, 'predicate')
            if ~strcmp(class(arg2{i}), 'predicate')
                result = false;
                return
            elseif ~comparepredicate(arg1{i}, arg2{i})
                result = false;
                return
            end
        else
            if ~comparevalue(arg1{i}, arg2{i})
                result = false;
                return
            end
        end
    end
    
    % if no mismatch is found, return true
    result = true;
end

function result = comparevalue( value1, value2 )
%cmpValue Compare 2 single values
    if isnumeric(value1) | islogical(value1)
        if ~isnumeric(value2) & ~islogical(value2)
            result = false;
            return
        else
            result = value1 == value2;
            return
        end
    else %ischar(pred1)
        if ~ischar(value2)
            result = false;
            return
        else
            result = strcmp(value1, value2);
            return
        end
    end
end

function result = eithernan( value1, value2 )
%eithernan Compare if either value is single NaN
    if length(value1) == 1 && ~iscell(value1) && isnan(value1)
        result = true;
        return
    elseif length(value2) == 1 && ~iscell(value2) && isnan(value2)
        result = true;
        return
    else 
        result = false;
        return
    end
end