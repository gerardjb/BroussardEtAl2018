function [Correction,x] = gammacorrection(Image,GammaValue)
% 
%Function to calculate the Gamma Correction for a particular image
%
%This function converts images or images files of monochromatic uint8
%class to an 8 bit, gamma corrected image.
% 

Err = 0;

%Determine if Image contains a filename or matrix
if ischar(Image)
    isfile = 1;
else
    isfile = 0;
end

%Determine if inputs look correct. Gamma default is 1.
if nargin < 2
    GammaValue = 1;
    disp('Default value for gamma = 1');
elseif nargin ==2 && GammaValue < 0
    GammaValue = 1;
    disp('GammaValue < 0, Default value considered, Gammavalue = 1');
elseif nargin > 2
    error('Error : Too many input parameters');
end

%Apply Gamma correction
if Err == 0 
    %Load image if Image is file, else declare Image as x
    if isfile
        x = imread(Image);
    else
        x = Image;
    end
    x = double(x);
    %Apply Gamma Correction to image
    Correction = (255 * (x/255).^ GammaValue);
    Correction = uint8(Correction);
    
end;