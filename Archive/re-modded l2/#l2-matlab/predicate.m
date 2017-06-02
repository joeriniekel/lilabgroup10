classdef predicate
    %PREDICATE Class to represent predicates of l2-models.
    %   p = predicate() creates an empty predicate/
    %   p = predicate('name{value1, 2}') or
    %   p = predicate('name', {'value1', 2}) create predicates with name
    %   'name' and arguments 'value1' and 2.
    %
    % PREDICATE Properties:
    %   name - The name of the predicate
    %   arg - Cell array with values for predicate arguments
    % 
    % See also ISPREDICATE, PREDCMP, PRED2STR, STR2PRED
    
    properties (SetAccess = 'private', GetAccess = 'public')
        name; % Predicate name
        arg;  % Cell array of arguments
    end
    
    methods
        function obj = predicate(varargin)
            switch nargin
                case 0 % empty model
                    obj.name =  '';
                    obj.arg =  {};
                case 1 % string input
                    obj = obj.string2predicate(varargin{1});
                case 2 % name and values input
                    obj.name =  varargin{1};
                    if iscell(varargin{2})
                        obj.arg =  varargin{2};
                    else
                        obj.arg =  varargin(2);
                    end
            end
        end
    end
    
    methods 
        % Conversion methods
        function result = double(a)
            try                   
                if isempty(a)
                    result = [];
                elseif numel(a) == 1 && length(a.arg) == 1
                    result = builtin('double', a.arg{1});
                elseif numel(a) > 1
                    [iMax, jMax] = size(a);
                    result = zeros(iMax,jMax);
                    for i=1:iMax
                        for j=1:jMax
                            result(i,j) = double(a(i,j));
                        end
                    end
                else
                    error('Conversion to double from predicate with multiple argumens is not possible.')
                end
            catch err
                throwAsCaller(err);
            end
        end
        
        function result = logical(a)
            try                   
                if isempty(a)
                    result = [];
                elseif numel(a) == 1 && length(a.arg) == 1
                    result = builtin('logical', a.arg{1});
                elseif numel(a) > 1
                    [iMax, jMax] = size(a);
                    result = zeros(iMax,jMax);
                    for i=1:iMax
                        for j=1:jMax
                            result(i,j) = logical(a(i,j));
                        end
                    end
                else
                    error('Conversion to logical from predicate with multiple arguments is not possible.')
                end
            catch err
                throwAsCaller(err);
            end
        end
        
        
        % Custom operators
        function result = plus(a,b)
            try                   
                result = double(a) + double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = minus(a,b)
            try                   
                result = double(a) - double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = uminus(a)
            try                   
                result = -double(a);
            catch err
                throwAsCaller(err);
            end
        end
        function result = uplus(a)
            try                   
                result = +double(a);
            catch err
                throwAsCaller(err);
            end
        end
        function result = times(a,b)
            try                   
                result = double(a) .* double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = mtimes(a,b)
            try                   
                result = double(a) * double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = rdivide(a,b)
            try                   
                result = double(a) ./ double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = ldivide(a,b)
            try                   
                result = double(a) .\ double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = mrdivide(a,b)
            try                   
                result = double(a) / double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = mldivide(a,b)
            try                   
                result = double(a) \ double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = power(a,b)
            try                   
                result = double(a) .^ double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = mpower(a,b)
            try                   
                result = double(a) ^ double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = lt(a,b)
            try                   
                result = double(a) < double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = gt(a,b)
            try                   
                result = double(a) > double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = le(a,b)
            try                   
                result = double(a) <= double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = ge(a,b)
            try                   
                result = double(a) >= double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = ne(a,b)
            try                   
                result = double(a) ~= double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = eq(a,b)
            try                   
                result = double(a) == double(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = and(a,b)
            try                   
                result = logical(a) & logical(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = or(a,b)
            try                   
                result = logical(a) | logical(b);
            catch err
                throwAsCaller(err);
            end
        end
        function result = not(a)
            try                   
                result = ~logical(a);
            catch err
                throwAsCaller(err);
            end
        end
    end
    
    methods (Hidden = true)
        function result = isempty(this)
            for i=1:numel(this)
                if ~isempty(this(i).name) || ~isempty(this(i).arg)
                    result = false;
                    return
                end
            end
            result = true;
        end
        
        function result = isnan(this)
            result = false;
        end      
        
        % Functions to convert string to predicate
         function result = string2predicate(this, string)
            % Convert a string of a predicate to Matlab variables
            valuesStart = min(strfind(string, '{'));
            if isempty(valuesStart)
                throw(MException('L2:STR2PRED:MissingBracket', ['Missing { for predicate ''' string '''']));
            elseif ~strcmp(string(end), '}')
                throw(MException('L2:STR2PRED:MissingBracket', ['Missing } for predicate ''' string '''']));
            else
                name = strtrim(string(1:valuesStart-1));
                if strcmp(name, '~')
                    name = NaN;
                end
                values = this.string2values(string(valuesStart+1:end-1));
                result = predicate(name, values);
            end
        end

        function result = string2values(this, string)
            % Split a comma separated list of values and convert each value
            try
                string = strtrim(string);
                result = {};

                while ~isempty(string)
                    [token, string] = strtok(string, ',');
                    token = strtrim(token);
                    string = strtrim(string);

                    % Token contains start of nested predicate
                    if strfind(token, '{')
                        closers = length(strfind(token, '{')) - length(strfind(token, '}'));
                        while closers > 0 & ~isempty(string)
                            [nextTok, string] = strtok(string, ',');
                            if strfind(nextTok, '{')
                                closers = closers+1;
                            end
                            if strfind(nextTok, '}')
                                closers = closers-1;
                            end
                            token = [token ',' nextTok];
                        end
                        result = [ result {this.string2predicate(token)} ];

                    % Otherwise, token contains a single construct    
                    else
                        result = [ result {this.string2value(token)} ];
                    end
                end
            catch err
                rethrow(err)
            end

        end

        function result = string2value(this, string)
        %singleValue Convert string to a single value
            try
                string = strtrim(string);

                [result,OK]=str2num(string);
                if OK && ~strcmp(string, 'i') && ~strcmp(string, 'j')
                    return
                elseif strcmp(string, '~')
                    result = NaN;
                else
                    result = string;
                end
            catch err
                rethrow(err)
            end
        end
        
       
    end
end

