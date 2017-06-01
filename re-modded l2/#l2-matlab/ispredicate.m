function tf = ispredicate( A )
%ISPREDICATE Determine whether item is a predicate array
%   This function returns logical 1 (true) if A is a predicate array and
%   logical 0 (false) otherwise.
%
%   See also PREDICATE, PREDCMP, PRED2STR, STR2PRED

    tf = isa(A, 'predicate');

end

