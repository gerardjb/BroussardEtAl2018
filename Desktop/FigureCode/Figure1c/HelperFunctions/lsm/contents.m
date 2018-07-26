% LSM File Toolbox.
% Version 1.1.0  12-Apr-2007
%
% Header files.
%   LSM_H  - Data structures used in reading in LSM info.
%   TIFF_H - Data structures used in reading in TIFF info.
%
% Public functions.
%   lsminfo      - Read in notes from Zeiss .lsm files.

% Private functions
%   channelinforead
%   ifdread             - Read IFD database from TIFF files according to TIFF standard.
%   scaninforead        - Read in the scaninfo database from Zeiss .lsm file
%   timestampsread
% 

%   Peter Li 12-Apr-07
%   Some rights reserved.  Licensed under Creative Commons: http://creativecommons.org/licenses/by-nc-sa/3.0/
