%Overlaying Red-Channel Derived Mask to Contrast-enhanced Images of Data
%by: Joey Broussard
%at UC Davis on 07/24/14
%
%This takes data from two files, one with a red mask (probably from mRuby)
%and green probe data, the other with the red FP data. It produces images
%of each dataset overlayed by a low opacity mask.

%% Clean Slate
CLEAN_SLATE = 1;

if CLEAN_SLATE
    clear all
    close all
end

%% Add Appropriate Libraries
importFunctions;

%% Toggles
USER_LEVEL = 1;%toggle to select isodata or threshold GUI

%% Open Files
%Open green data with red mask file, Package it
[file1,path1] = uigetfile('*.*',...
    'Select the file containing the Green Data and Mask');
fname1 = fullfile(path1,file1);
[mPix,nPix,nImage,Img1] = openTiff(fname1,2);
Gdata1 = Img1.green;
Rdata1 = Img1.red;

%Open file with red data It, Package it
[file2,path2] = uigetfile('*.*',...
    'Select the file with the Red Data');
fname2 = fullfile(path2,file2);
[mPix,nPix,nImage,Img2] = openTiff(fname2,2);
Rdata2 = Img2.red;

%% Create Axon and Soma Masks
%Make maxintensity projection for red mask data
redMax1 = max(Rdata1,[],3);

%Mask this data using isodata algorithm or user-defined
if USER_LEVEL
    redMaskPrep = threshPick(redMax1);
else
    [level,redMaskPrep] = isodata(redMax1);
end

uiwait(gcf);
%Clean up small ROIs
smallROISize = 100;
redMaskPrep2 = removeSmallROI(redMaskPrep,smallROISize);


%User Select ROI containing the cell

imshow(redMax1)
fontSize = 10;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
message = sprintf(...
    'Make an ROI around the axon.\nDouble-click inside polygon to accept it.');
% Display instructions as title.
title(message, 'FontSize', fontSize);
%uiwait(msgbox(message));
userMask = roipoly();
close gcf;

%Combine User and Isodata collected masks
redMask = redMaskPrep2&userMask;

%% Contrast Enhance Images
%Make max intensity projections for both red and green data
greenMax1 = max(Gdata1,[],3);
redMax2 = max(Rdata2,[],3);

%Apply gamma correction to both images
GammaGreen = 0.6;
GEnhance = gammacorrection(greenMax1,GammaGreen);
GammaRed = 0.6;
REnhance = gammacorrection(redMax2,GammaRed);

%Change Colors from gray to respective image color
GEnhance_color = gray2oneColor(GEnhance,'g');
REnhance_color = gray2oneColor(REnhance,'r');

%% Make Overlay Figures

%Figure for mask creation
figure(1)
subplot(1,2,1)
imshow(redMax1)
title('Data used to Create Mask')

subplot(1,2,2)
imshow(redMask)
title('Mask created from mRuby High Exposure')

%Figure with Green Data and mask overlay
figure(2)
subplot(1,2,1)
imshow(greenMax1);
title('Maximum intensity Projection of Green Data')

subplot(1,2,2)
imshow(GEnhance_color)
showMaskAsOverlay2(0.1,redMask,'w');
titleStrG = sprintf('Mask overlayed to the Green Data with Gamma Value %0.2g'...
    ,GammaGreen);
title(titleStrG);

%Figure with Red Data and mask overlay
figure(3)
subplot(1,2,1)
imshow(redMax2);
title('Maximum intensity Projection of Red Data')

subplot(1,2,2)
imshow(REnhance_color)
showMaskAsOverlay2(0.1,redMask,'w',[],0);
titleStrR = sprintf('Mask overlayed to the Red Data with Gamma Value %0.2g'...
    ,GammaRed);
title(titleStrR);