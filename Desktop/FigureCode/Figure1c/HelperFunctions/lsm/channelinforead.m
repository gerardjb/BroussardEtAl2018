function lsminf = channelinforead(fid, lsminf)
%CHANNELINFOREAD Read in channelinfo from Zeiss LSM file format
%   LSMINF = CHANNELINFOREAD(FID, LSMINF)
%
%   FID is the fopen file identifier to read data from.  LSMINF is the
%   previously read LSM headers up to this point, including file offset
%   locations for various channel info data.
%

%   Peter Li 11-Apr-07
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/

error(nargchk(2, 2, nargin));

LSM_H; % CHANNELCOLORS from header file

if lsminf.OffsetChannelColors == 0,
    lsminf(1).ChannelColors = [];
    return
else
    if fseek(fid, lsminf.OffsetChannelColors, 'bof'),
        error(['Received error on file seek to OffsetChannelColors(' lsminf.OffsetChannelColors '): ' ferror(fid)]); 
	end
	
	channelcolors = structread(fid, CHANNELCOLORS);
	
	% Read in Channel RGB values
	fseek(fid, lsminf.OffsetChannelColors + channelcolors.ColorsOffset, 'bof');
	for i = 1:channelcolors.NumberColors
        R = fread(fid, 1, 'uint8');
        G = fread(fid, 1, 'uint8');
        B = fread(fid, 1, 'uint8');
        channelcolors.Colors{i} = [R G B];
	end
	
	% Read in Channel Names
	fseek(fid, lsminf.OffsetChannelColors + channelcolors.NamesOffset, 'bof');
	for i = 1:channelcolors.NumberNames
        namelength = fread(fid, 1, 'uint32');
        name = char(fread(fid, namelength, 'char')');
        if uint8(name(end)) == 0
            name = name(1:end-1);
        end
        channelcolors.Names{i} = name;
	end
	
	lsminf(1).ChannelColors = channelcolors;
end
