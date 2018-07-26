function scaninfo = scaninforead(fid, offset)
%SCANINFOREAD Read in scaninfo database of Zeiss LSM file format
%   SCANINFO = SCANINFOREAD(FID, [OFFSET])
%
%   OFFSET should be the file position for the start of the scaninfo database.
%   If FID is already at the correct file position, OFFSET should be omitted.
%

%   Peter Li 30-Aug-05
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(1, 2, nargin));

LSM_H; % SCANINFO_HEXTAGMAP from header file
scaninfo = struct([]);

% Seek to the scaninfo database file position if an offset is given, otherwise
% leave FID in its current position.
if nargin > 1
    if fseek(fid, offset, 'bof'), 
        error(['Received error on file seek to SCANINFO_OFFSET(' offset '): ' ferror(fid)]); 
    end
end

% The algorithm for reading the scaninfo database depends on keeping track of a
% "level" hierarchy.  As the database is read, some entries are level
% instructions, navigating up or down the hierarchy.  The database is done when
% the hierarchy steps back to level 0.
level = 0;
while (1)
    taghex      = dec2hex(fread(fid, 1, 'uint32'));
    typecode    = fread(fid, 1, 'uint32');
    size        = fread(fid, 1, 'uint32');
    tag = structmap(['h' taghex], SCANINFO_HEXTAGMAP);

    switch typecode
        case 0 % Special case: this is a level instruction entry
            if (hex2dec(taghex) == hex2dec('FFFFFFFF'))
                level = level - 1;
            else
                level = level + 1;
            end
        case 2 % string
            count = size;
            value = char(fread(fid, count, 'uchar')');
            value = value(1:end-2);
        case 4 % int32
            count = size / 4;
            value = fread(fid, count, 'uint32');
        case 5 % float64
            count = size / 8;
            value = fread(fid, count, 'float64');
        otherwise
            fseek(fid, size, 0);
            value = '?';
    end

    % If this was just a level instruction entry ignore it, otherwise try to
    % record entry.
    if typecode > 0
        if ~isempty(tag)
            scaninfo(1).(tag) = structappend(scaninfo, tag, value);
        else
            scaninfo(1).(['h' taghex]) = structappend(scaninfo, ['h' taghex], value);
        end
    end

    if (level == 0)
        break;
    end
end
