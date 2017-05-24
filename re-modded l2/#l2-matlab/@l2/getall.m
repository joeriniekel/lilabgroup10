function result = getall( trace, t, varargin )
%GETALL Get a particular set of values for a predicate. 
%   The first two arguments contain the trace and timepoint t to be
%   checked. There are three ways of using this function.
%   1) The third argument is a predicate containing the target pattern
%   2) The third argument is a string pattern.
%   3) The third argument is the targeted predicate and the fourth the
%   pattern of values only.
%
%   Example:
%   l2.GETALL( trace, 1, 1x1 predicate)
%   and is similar to
%   l2.GETALL( trace, 1, 'pred{value1, ~}')
%   and is similar to
%   l2.GETALL( trace, 1, 'pred', '{value1, ~}')
%
%   See also EXISTS, PREDICATE

    % Retrieve predicate and pattern, depending on usage of function
    if nargin==3
        if ischar(varargin{1})
            pattern = predicate(varargin{1});
        else
            pattern = varargin{1};
        end
    else
        pattern = predicate(varargin{1}, varargin{2});
    end
    
    % Get all values for that predicate
    pAll = [trace(t).(pattern.name)];
    
    % if no values exist for predicate, return empty result
    if isempty(pAll)
        result = pAll;
        return
    end
    
    result = [];
    for i=1:length(pAll)
        if predcmp(pAll(i), pattern)
            result = [result pAll(i)];
        end
    end
end