function newvalue = structappend(struct, field, value)
%STRUCTAPPEND
%   VALUE = STRUCTAPPEND(STRUCT, FIELD, VALUE)
%
%   When adding data to struct, gracefully handle duplicated field names by
%   concatenating values into a cell array.
%

%   Peter Li 30-Aug-05
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(3, 3, nargin));

if isfield(struct, field) % This field already initialized!
    if iscell(struct.(field)) % Already converted this field into cells
        newvalue = struct.(field);
        newvalue{end+1} = value;
    else                      % Otherwise, convert this field into cells
        newvalue = {struct.(field), value};
    end
else % Make a new field 
    newvalue = value;
end
