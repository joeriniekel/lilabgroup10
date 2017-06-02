function [exists, values] = exists( trace, t, varargin )
%EXISTS Check whether a particular value for a predicate exists. 
%   The first two arguments contain the trace and timepoint t to be
%   checked. There are two ways of using this function.
%   1) The third argument is the targeted predicate and lastly a 
%   desired predicate pattern is provided. 
%   2) The third argument is a string pattern, including the predicate.
%   Returned is primarly a boolean indicating whether or not the targeted 
%   pattern is present, but the matching values are returned as well.
%
%   Example:
%   [exists, values] = l2.EXISTS( trace, 1, 'pred', {'value1', NaN})
%   and is similar to
%   [exists, values] = l2.EXISTS( trace, 1, 'pred{value1, NaN}')
%
%   See also GETALL

    % Use getall to find matching value and derive exists from that.
    values = l2.getall(trace, t, varargin{:});
    exists = ~isempty(values);
end
