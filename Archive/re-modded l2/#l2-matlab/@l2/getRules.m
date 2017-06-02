function fncs = getRules()
%GETRULES Return list of local function handles
%   Detailed explanation goes here
    fncs = {};
    %mfilename('fullpath')
    ST = dbstack('-completenames');
    filename = ST(2).file;
    fid = fopen(filename);
    [~, nRules, ~]=fileparts(filename);
    
    tline = fgetl(fid);
    while ischar(tline)
        tline = strtrim(tline);
        
        %remove comments (if block comment, removes entire block)
        tline = removeComments(tline, fid);
        
        %if line is empty, skip
        if isempty(tline)
            tline = fgetl(fid);
            continue
        else %start of a function
            [fName tline] = processFunction(tline, fid);
           if ~strcmp(fName, nRules)
               fncs{ length(fncs)+1 } = fName;
           end
        end
    end
    fclose(fid);
end

function tline = removeComments(tline, fid)

    %if line is empty, done!
    if isempty(tline)
        return
    end
    
    %if blockcomment, skip it (always on empty line)    
    if length(tline) == 2 && strcmp(tline,'%{')
        skip = 1;
        while skip > 0
            tline = strtrim(fgetl(fid));
            if length(tline) == 2 && strcmp(tline,'%}')
                skip = skip - 1;
            elseif length(tline) == 2 && strcmp(tline,'%{')
                skip = skip + 1;
            end
        end
        tline = '';
        return;
    end
    
    possibleComments = strfind(tline, '%');
    for index = possibleComments
        % if enclosed by quotes, ignore this sign
        quotes = strfind(tline, '''');
        if ~isempty(quotes)
            % if uneven amount of quotes left of '%', then enclosed
            if mod(sum(quotes<index),2)
                continue
            end
        end 
        
        if length(tline) > index    % more characters after '%'
            %remove rest of tline and stop
            tline = strtrim( tline(1:index-1) );
            break
        else                        % '%' is last sign
            if index == 1           % '%' is only sign on the line
                tline = '';
                return
            else                    % return tline up to index
                tline = strtrim( tline(1:index-1) );
                return
            end
        end
    end

end


function [fName tline] = processFunction(tline, fid)
    % get function name, tline starts with:
    % function = <name>(<arguments>)
    
    leftIndex = strfind(tline, '=');
    rightIndex = strfind(tline, '(');
    endIndex = strfind(tline, ')');
    fName = strtrim( tline(leftIndex(1)+1:rightIndex(1)-1) );
    
    if length(tline) > endIndex
        tline = tline(endIndex:end);
    else
        tline = fgetl(fid);
    end
    
    %find end of function, returns remainder of that line
    tline = findEnd(tline, fid);
end
        
function tline = findEnd(tline, fid)
    while ischar(tline)
        tline = strtrim(tline);
        
        %remove comments (if block comment, removes entire block)
        tline = removeComments(tline, fid);
        
        %if line is empty, skip
        if isempty(tline)
            tline = fgetl(fid);
            continue
        else %check for end or nested arguments that require end
            args = {'function' 'for' 'while' 'switch' 'try' 'if' 'parfor'};
            match = cellfun(@(x) regexp(tline,['\<' x '\>']), args, 'UniformOutput', false);
            possibleNested = [match{:}];
            possibleEnds =  regexp(tline,'\<end\>');
           
            %check if not between quotes
            quotes = strfind(tline, '''');
            if ~isempty(quotes)
                nested = [];
                ends = [];
                for index = possibleNested
                    % if even amount of quotes left of '%', then not enclosed
                    if ~mod(sum(quotes<index),2)
                        nested = [nested index];
                    end
                end
                for index = possibleEnds
                    % if even amount of quotes left of '%', then not enclosed
                    if ~mod(sum(quotes<index),2)
                        ends = [ends index];
                    end
                end
            else
                nested = possibleNested;
                ends = possibleEnds;
            end
            
            % check if end is not between brackets
            brackets = [strfind(tline, '(') strfind(tline, ')')];
            if ~isempty(brackets)
                possibleEnds = ends;
                ends = [];
                for index = possibleEnds
                    % if even amount of quotes left of '%', then not enclosed
                    if ~mod(sum(brackets<index),2)
                        ends = [ends index];
                    end
                end
            end
            
            % if there are no ends or nested arguments, continue on next line
            if isempty(nested) && isempty(ends)
                tline = fgetl(fid);
                continue
            % if line contains both, see which one comes first and forget about the others
            elseif ~isempty(nested) && ~isempty(ends)
                if min(nested) < min(ends)
                    ends = [];
                else 
                    nested = [];
                end
            end
                
                
            % if line only contains ends, remove first one and return rest of tline
            if ~isempty(ends)
                if ends(1)+3 < length(tline) % first end is not end of line
                    tline = strtrim( tline(ends(1)+3:end) );
                    return
                else                    % return empty line
                    tline = '';
                    return
                end       
                
            % if line contains nested, remove first one and search for its end    
            else
                %first remove stuff before the argument
                if nested(1)~=1
                    tline = tline(nested(1):end);
                end
                
                %remove the argument itself
                argEnd = strfind(tline, ' ');
                if isempty(argEnd)
                    tline = fgetl(fid);
                    tline = findEnd(tline, fid);
                    continue
                else
                    tline = tline(argEnd(1):end);
                    tline = findEnd(tline, fid);
                    continue
                end
            end
        end
    end
end

%{

           %find end of function
           ends = 1;
           tline = fgetl(fid);
           while ends > 0 && ischar(tline)
               tline = strtrim(tline);


               %if line is comment or block comment
               if strcmp(tline(1),'%')
                   %if length tline == 1 or not block comment
                   if length(tline) == 1 || ~strcmp(tline(1:2),'%{')
                       tline = fgetl(fid);
                       continue
                   else
                       skip = true;
                       while skip
                            tline = strtrim(tline);
                            skip = ~strcmp(tline(1:2),'%}');
                            tline = fgetl(fid);
                       end
                       continue
                   end
               end
               
  

               %if line contains argument requiring end, ends + 1
               args = {'function' 'for' 'while' 'switch' 'try' 'if' 'parfor'};
               match = cellfun(@(x) regexp(tline,['\<' x '\>']), args, 'UniformOutput', false);
               if ~isempty([match{:}])
                   ends = ends + 1;
               end

               %if line contains an end, ends - 1
               if regexp(tline,'\<end\>')
                   ends = ends -1;
               end

               %go to next line
               tline = fgetl(fid);
           end
           
        %block comment outside a function
        elseif length(tline) >= 2 && strcmp(tline(1:2),'%{') skip = true;
           skip = true;
           while skip
                tline = strtrim(tline);
                skip = length(tline) < 2 || ~strcmp(tline(end-1:end),'%}');
                tline = fgetl(fid);
           end
           continue
        else
            tline = fgetl(fid);
        end
    end
    fclose(fid);

end
%}