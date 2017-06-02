function result = str2pred( string )
%STR2PRED Convert string in predicate format
%   Convert a given string to predicate format. Input can only be a string
%   containing a single predicate. If conversion fails, an empty predicate
%   is returned.
%
%   See also PREDICATE, ISPREDICATE, PREDCMP, PRED2STR

    if ~ischar(string)
        throwAsCaller(MException('L2:STR2PRED:unrecognizedFormat', ['Unable to convert ' class(string) ' to predicate, use string.']));
    else
        string = strtrim( string );
    end
    
    try
        result = predicate(string);
    catch
        result = predicate();
    end
end