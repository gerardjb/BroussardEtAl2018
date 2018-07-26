function readstruct = structread(fid, structdef);
% STRUCTREAD read in structured data
%   READSTRUCT = STRUCTREAD(FID, STRUCTDEF)
%
%   Read in structured data from file handle FID according to STRUCTDEF,
%   which is in the ordered format of structdef.fieldname = datatype
%   For example, structdef.isbn = '10*uint8' means that the ISBN field to be
%   read consists of 10 unsigned integers.
%

%   Peter Li 30-Aug-05
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(2, 2, nargin));
if ~isstruct(structdef), error('Second argument must be a struct definition'); end

readstruct = struct([]);

% Get file length
pos = ftell(fid);
fseek(fid, 0,'eof');
flen = ftell(fid);
if fseek(fid, pos, 'bof') == -1, error(['Received error on file seek: ' ferror(fid)]); end

% Calculate how many bytes the file has left, compare to size of struct
fleft = flen - pos;
structsize = structsizeof(structdef);
if structsize > fleft, error(['Size of Struct (' num2str(structsize) ') > Bytes Left in File (' num2str(fleft) ')']); end

% Main section   
fields = fieldnames(structdef);
for i = 1:length(fields),
    field = fields{i};
    datatype  = structdef.(field);

    % Parse datatype
    starind = strfind(datatype, '*');
    if length(starind) > 1, error(['Incorrect datatype format: ' datatype]); end
    if starind,
        num = str2num(datatype(1:starind-1));
        datatype = datatype(starind + 1:end);
    else
        num = 1;
    end

    % Read value(s)
    [value, readnum] = fread(fid, num, datatype);
    if readnum ~= num, error(['Failed to read more than ' num2str(readnum) ' values for ' field '(' num2str(num) ')']); end

    % If value is character string, convert to char.
    if strfind(datatype, 'char'),
        value = char(value);
    end

    % Handle fieldnames ending in numbers, e.g. if fields UNKNOWN1, UNKNOWN2,
    % etc. exist, we will concatenate them all as a cell array into a single
    % UNKNOWN field.

    % See if filename ends in a number
    for i = fliplr(1:length(field)),
        if ~length(str2num(field(i:end))),
            break;
        end
    end
    if i < length(field),
        appendnum = str2num(field(i+1:end));
        appendfield = field(1:i);

        % See if appendfield exists and whether appendnum matches the number
        % of elements to be stored there.
        if ~isfield(readstruct, appendfield),
            if appendnum == 1, field = appendfield; end
        else
            if ~iscell(readstruct.(appendfield)),
                if appendnum == 2, field = appendfield; end
            else
                if appendnum == length(readstruct.(appendfield)) + 1, field = appendfield; end
            end
        end
    end

    readstruct(1).(field) = structappend(readstruct, field, value');
end
