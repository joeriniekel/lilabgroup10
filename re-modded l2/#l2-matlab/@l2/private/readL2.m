function structs = readL2( filename, file )
%READL2 Read a file in the l2 format
%   Read a file with the extension .l2. The file consists one or more sets,
%   where each set contains a one or more lines with two values. If no set
%   is specified, the contents are returned as a set named 'default'.
%
%   <setname> (
%       <name>; <value(s)>
%       <name>; <value(s)>
%   )

    try
        fid = fopen(filename);
    
        %read using both the set delimeter () and name-value delimiter ;
        %set name is in first column, where second column is empty
        C = textscan(fid,'%s %s','Delimiter',';()', 'CommentStyle', '%', 'CollectOutput', true);
        contents = C{1};

        structs = struct;
        setName = '';
        values = struct;
        for i=1:size(contents,1)
            content1 = contents{i,1};
            content2 = contents{i,2};
            
            %if first column is empty, second column must be empty as well
            if isempty(content1)
                if ~isempty(content2)
                    throw(MException('L2:READL2:UnexpectedText', 'Unexpected text ''%s'' in ''%s''.', content2, filename));
                else
                    continue
                end
            else
                content1 = strtrim(content1);
                if ~isempty(content2)
                    content2 = strtrim(content2);
                end
            end

            % If setName is still empty, and first line contains a
            % name-value pair, initialise setName as default
            if isempty(setName) & ~isempty(content1) & ~isempty(content2)
                setName = 'default';
            end
            
            %if only second column is empty, start of a newset
            if ~isempty(content1) && isempty(content2)
                %save results (if any)
                if (~isempty(fieldnames(values)))
                    try
                        structs.(setName) = values;
                    catch err
                        throwAsCaller(MException('L2:READL2:InvalidFieldName', 'In ''%s'' invalid set: ''%s''.\n\tNo hyphens (-), slashes (/) or white space ( ) allowed.', filename, setName));
                    end
                end
                %change setName and clear values
                setName = content1;
                values = struct;

            %line contains name value pair    
            else
                
                if strcmp(file, 'scenarios')
                    try
                        pred = predicate(content1);
                        if isfield(values, pred.name)
                            values.(pred.name) = {values.(pred.name){:} {str2other(content2), pred}};
                        else
                            values.(pred.name) = {{str2other(content2), pred}};
                        end
                    catch err
                        if exist('pred','var')
                            throwAsCaller(MException('L2:READL2:InvalidFieldName', 'In ''%s'' invalid predicate name: ''%s''.\n\tNo hyphens (-), slashes (/) or white space ( ) allowed.', filename, pred.name));
                        else
                            throwAsCaller(MException('L2:READL2:UnexpectedText', 'In ''%s'' invalid predicate ''%s'': \n\t%s.', filename, content1, err.message));
                        end
                    end
                else
                    try    
                        values.(content1) = str2other(content2);
                    catch err
                        throwAsCaller(MException('L2:READL2:InvalidFieldName', 'In ''%s'' invalid name: ''%s''.\n\tNo hyphens (-), slashes (/) or white space ( ) allowed.', filename, content1));
                    end
                        
                end
            end
        end

        %save the last set
        try
            structs.(setName) = values;
        catch err
            throwAsCaller(MException('L2:READL2:InvalidFieldName', 'In ''%s'' invalid set: ''%s''.\n\tNo hyphens (-), slashes (/) or white space ( ) allowed.', filename, setName));
        end
            
    catch err
        % if given file is not found
        if strcmp(err.identifier, 'MATLAB:FileIO:InvalidFid')
            throwAsCaller(MException('L2:READL2:FileNotFound', 'Unable to find ''%s''.', filename));
        
        % if reading the file, unexpected text shows up
        elseif strcmp(err.identifier, 'L2:READL2:UnexpectedText')
            throwAsCaller(err);
            
        % if an invalid set name is provided
        elseif strcmp(err.identifier, 'L2:READL2:InvalidFieldName')
            throwAsCaller(err)
            %throwAsCaller(MException('L2:READL2:InvalidFieldName', 'In ''%s'' invalid field name: ''%s''.\n\tNo hyphens (-), slashes (/) or white space ( ) allowed.', filename, content1));
        
        else
            err.identifier
            rethrow(err)
        end
    end
end

function result = str2other(value)
%STR2OTHER convert a value to another type
    value = strtrim(value);
    
    % if it's empty, return value as is
    if isempty(value)
        result = value;
        return
    end
    
    % check if it's a predicate
    %pred = str2pred(value);
    %if ~isempty(pred)
    %    result = pred;
    %    return
    %end     
    
    % check if it's a cell array
    if strcmp(value(1), '{') && strcmp(value(end), '}')
        result = strsplit(value(2:end-1), ',');
        for i=1:length(result)
            result{i} = str2other(result{i});
        end
        return
    end
    
    % if value is a single numeric value
    [result,OK]=str2num(value);
    if OK && ~strcmp(value, 'i') && ~strcmp(value, 'j')
        return
    end
    
    % if value is a single string (no comma's)
    result = value;
end