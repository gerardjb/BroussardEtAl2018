function [ifd,offset] = ifdread(fid)
% IFDREAD read IFD entries from TIFF file according to TIFF standard
%   [IFD,OFFSET] = IFDREAD(FID)
%
%   The idf database is read into struct IFD.  The OFFSET to the beginning of
%   image data is the last  piece of information in the database.
%
%   Peter Li 30-Aug-05
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(1, 1, nargin));

TIFF_H; % ENTRY_LENGTH, IFD_TAGMAP from header file

if fseek(fid, 8, 'bof'), error(['Received error on file seek to offset 8: ' ferror(fid)]); end
num_entries = fread(fid, 1, 'uint16');
entries_start_pos = ftell(fid);

for i = 1:num_entries,
    % Go to ifd entry i
    if fseek(fid, entries_start_pos + ((i - 1) * ENTRY_LENGTH), 'bof'), error(['Received error on file seek to position for entry ' num2str(i) ': ' ferror(fid)]); end

    ifd(i).tagcode     = fread(fid, 1, 'uint16');
    ifd(i).typecode    = fread(fid, 1, 'uint16');
    ifd(i).count       = fread(fid, 1, 'uint32');
    ifd(i).tag = structmap(ifd(i).tagcode, IFD_TAGMAP);
    switch ifd(i).typecode
        case {1, 2, 6, 7}   % uint8, uchar, int8, undef8
            ifd(i).bytes = 1;
        case {3, 8}         % uint16, int16
            ifd(i).bytes = 2;
        case {4, 9, 11}     % uint32, int32, float32
            ifd(i).bytes = 4;
        case {5, 10, 12}    % uint32 / uint32, int32 / int32, float64
            ifd(i).bytes = 8;
        otherwise
            error(['Illegal typecode: ' num2str(ifd(i).typecode) '.  Misformed IFD']);
    end

    % If the data take up more than 4 bytes, then the value is a
    % pointer to the data.  Otherwise the value is the data.
    if ((ifd(i).bytes * ifd(i).count) > 4)
        ifd(i).offset = fread(fid, 1, 'uint32');
        if fseek(fid, ifd(i).offset, 'bof')
            error(['Received error on file seek to data for entry ' num2str(i) '(' num2str(ifd(i).offset) '): ' ferror(fid)]);
        end
    end

    switch ifd(i).typecode
        case 1 % uint8
            ifd(i).value = fread(fid, ifd(i).count, 'uint8');
        case 2 % uchar
            ifd(i).value = fread(fid, ifd(i).count, 'uchar');
            char(ifd(i).value);
        case 3 % uint16
            ifd(i).value = fread(fid, ifd(i).count, 'uint16');
        case 4 % uint32
            ifd(i).value = fread(fid, ifd(i).count, 'uint32');
        case 5 % Two uint32, first a numerator, then a denominator, representing a fraction
            ifd(i).value = fread(fid, 2 * ifd(i).count, 'uint32');
            numerators   = fid(i).value(1:2:end);
            denominators = fid(i).value(2:2:end);
            ifd(i).value = double(numerators) ./ double(denominators);    
        case 6 % int8
            ifd(i).value = fread(fid, ifd(i).count, 'int8');
        case 7 % undef8
            ifd(i).value = fread(fid, ifd(i).count, 'uint8');
        case 8 % int16
            ifd(i).value = fread(fid, ifd(i).count, 'int16');
        case 9 % int32
            ifd(i).value = fread(fid, ifd(i).count, 'int32');
        case 10 % Two int32, first a numerator, then a denominator, representing a fraction
            ifd(i).value = fread(fid, 2 * ifd(i).count, 'int32');
            numerators   = fid(i).value(1:2:end);
            denominators = fid(i).value(2:2:end);
            ifd(i).value = double(numerators) ./ double(denominators);    
        case 11 % float32
            ifd(i).value = fread(fid, ifd(i).count, 'float32');
        case 12 % float64
            ifd(i).value = fread(fid, ifd(i).count, 'float64');
        otherwise
            error(['Illegal typecode: ' num2str(ifd(i).typecode) '.  Misformed IFD']);
    end
end

fseek(fid, entries_start_pos + (num_entries * ENTRY_LENGTH), -1);
offset = fread(fid, 1, 'uint32');
