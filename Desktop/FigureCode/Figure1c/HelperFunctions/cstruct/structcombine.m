function combined = structcombine(struct1, struct2)
%STRUCTCOMBINE
%   COMBINED = STRUCTCOMBINE(STRUCT1, STRUCT2)
%
%   Crash structs together, using STRUCTAPPEND to resolve overlapping fields.
%   This will take all the fieldnames of STRUCT2 and add them to STRUCT1, then
%   copy all the data of STRUCT2 into the new expanded STRUCT1, including data
%   in higher dimensions of struct2.
%
%   For example, if STRUCT1 is:
%       struct1.x = 0
%
%   And STRUCT2 is:
%       struct2.y = 0
%       struct2(50,50).y = 0
%
%   Then the result, COMBINED is a 50-by-50 struct array with fields x and y and
%   almost no data except for fields combined(1).x, combined(1).y, and
%   combined(50,50).y, all of which are set to zero.
%
%   For example of overlapping fields, if:
%       struct1.x = 0
%       struct2.x = 0
%   
%   Then combined.x will be cell array {0 0}.

%   Peter Li 12-Apr-07
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(2, 2, nargin));

% Start with just STRUCT1
combined = struct1;

% Iterate through the fields of STRUCT2, adding them to STRUCT1
fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);
nlength2 = prod(size(struct2));
for i = 1:length(fields2)
    % Check that the field name we're about to add doesn't already exist in
    % STRUCT1
    if ~length(strmatch(fields2{i}, fields1))
        [combined(1:nlength2).(fields2{i})] = deal(struct2.(fields2{i}));
    else
        combined(nlength2).(fields2{i}) = [];
        for j = 1:nlength2
            combined(j).(fields2{i}) = structappend(combined(j), fields2{i}, struct2(j).(fields2{i}));
        end
    end
end
