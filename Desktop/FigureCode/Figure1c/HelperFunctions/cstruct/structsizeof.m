function [size, errormsg] = structsizeof(structdef)
%STRUCTSIZEOF Determine size (in bytes) of struct to be read
%   SIZE = STRUCTSIZEOF(STRUCTDEF)
%

%   Peter Li 30-Aug-05
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

if nargin ~= 1 | ~isstruct(structdef), error('First argument must be struct definition'); end

size = 0;
errormsg = '';
fields = fieldnames(structdef);
for i = 1:length(fields);
    datatype = structdef.(fields{i});

    % Parse num
    starind = strfind(datatype, '*');
    if length(starind) > 1, error(['Incorrect datatype format: ' datatype]); end
    if starind
        num = str2num(datatype(1:starind(1) - 1));
        datatype = datatype(starind(1)+1:end);       
    else
        num = 1;
    end
    
    switch datatype
        case {'uchar' 'schar' 'uint8' 'int8'}
            size = size + num;
        case {'uint16' 'int16'}
            size = size + (2 * num);
        case {'uint32' 'int32' 'single' 'float32'}
            size = size + (4 * num);
        case {'uint64' 'int64' 'double' 'float64'}
            size = size + (8 * num);
        otherwise
            size = -1;
            errormsg = ['Unknown datatype: ' datatype];
            return;
    end
end
