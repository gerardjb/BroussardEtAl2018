function [out,ADRout] = CalcADR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating ADR and Overlaying Red-Channel Derived Mask to Images of Data
%Joey Broussard
%Tian Lab, UC Davis
%09/08/14
%
%This takes data from two files, one with a red mask (probably from mRuby)
%and green probe data, the other with the red FP data. It produces images
%of each dataset overlayed by a low opacity mask. It also produces an image
%of the cell with pixels binned as a function of distance to the
%somatodendritic compartment.
%
%Currently, it is set to analyze the attached sample file in the directory
%where the file is located. However, if the DEBUG variable is set to 0, the
%user will be promepted to select a file to analyze.
%
%   Outputs:
%       out - double with unnormalized ADR, normalized ADR, ADR for the red
%           fluorophore, laser power for the 488 and 561 lines, and the
%           size of the axonal arbor in micron^2
%       ADRout - structure containing cell identifiers, and data used to
%           plot the different ADR metrics as a function of distance from
%           the somatodenritic comparment
%
%   Plots:
%       1 - Colormap of cell indicating where spatial bins located
%       2 - Scatterplot of uADR as a function of distance from soma with
%           same colormap as in 1. Also shows number of pixels in each
%           spatial bin.
%       3 - Scatterplot of nADR as a function of distance from soma with
%           same colormap as in 1. Also shows number of pixels in each
%           spatial bin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Open Files, load data
clear all
close all

%Pull in help functions
path1 = fileparts(which('calcADR.m'));
addpath(genpath([path1,'/HelperFunctions']));
USER_LEVEL = 1;%toggle user defined threshold or iso mask, 1=user
DEBUG = 1;%Reuse a selected file to debug

%Open green data with red data, Package it
if DEBUG
    file1 = '20xG&R2.tif';
    fname1 = fullfile(path1,file1);
else
    [file1,path1] = uigetfile('*.*',...
        'Select the file containing the Green Data and Red Data');
    fname1 = fullfile(path1,file1);
end
Img1 = openTiff(fname1,2);
Gdata1 = Img1.green;
Rdata2 = Img1.red;

%Open file with mask, Package it
file2 = [file1(1:end-4),'mask.tif'];
fname2 = fullfile(path1,file2);
Img2 = openTiff(fname2,2);
Rdata1 = Img2.red;

%Get laser powers for this image set
lsmfname1 = strgsub(fname1,'tif','lsm');
[lsm1,~,~] = lsminfo(lsmfname1);
lasers = lsm1.ScanInfo.POWER;
if iscell(lasers)
    laserPower488 = lasers{2};%green laser is in the second position
    laserPower561 = lasers{1};
else
    laserPower488 = lasers;%green laser only one used
    laserPower561 = [];%red laser unecessary due to high red fluorophore expression
end

%Get scale of the pixels for this image to calculate axon pixels binned
%to different distances from the soma.
targetRadius = 50;%radius of each annulus in microns
scale = lsm1.VoxelSizeX*10^6;%pixel edge in microns
radius = targetRadius/scale;%radius for concentric circles in pixels


%% Create unique identifier series for each processed file
%Pull data about this experiment from the path
splitPath = strsplit(path1,'/');plate = splitPath{end - 1};
splitStub = strsplit(file1,'.');cell = splitStub{1};

%Populate structure ID fields
%Date - Note: removed because directories do not specify date
%id_date = splitPath{end - 3};
%newADRout.id_date = id_date;
%Construct
id_const = splitPath{end - 2};
newADRout.id_const = id_const;
%CellID
id_cell = sprintf('P%sC%s',plate(end),cell(end));
newADRout.id_cell = id_cell;

%% Create Axon and Soma Masks
%Make maxintensity projection for red mask data
redMax1 = max(Rdata1,[],3);

%Mask this data using isodata algorithm or user-defined
if USER_LEVEL
    redMaskPrep = threshPick(redMax1);
else
    [level,redMaskPrep] = isodata(redMax1);
end

%Clean up small ROIs
smallROISize = 200;
redMaskPrep = removeSmallROI(redMaskPrep,smallROISize,8);


%User Select ROI containing the whole cell
imshow(redMax1)
fontSize = 10;
set(gcf, 'units','normalized','outerposition',[0 0 1 0.9]);
message = sprintf(...
    'Make an ROI around the whole cell.\nDouble-click inside polygon to accept it.');
% Display instructions as title.
title(message, 'FontSize', fontSize);
%uiwait(msgbox(message));
cellMaskPrep = roipoly();
close gcf;

%User Select ROI containing the soma
imshow(redMax1)
set(gcf, 'units','normalized','outerposition',[0 0 1 0.9]);
message = sprintf(...
    'Make an ROI around the soma.\nDouble-click inside polygon to accept it.');
% Display instructions as title.
title(message, 'FontSize', fontSize);
%uiwait(msgbox(message));
somaMaskPrep = roipoly();
close gcf;

%User Select ROI containing background pixels, get background value
imshow(redMax1)
set(gcf, 'units','normalized','outerposition',[0 0 1 0.9]);
message = sprintf(...
    'Make an ROI around the background.\nDouble-click inside polygon to accept it.');
% Display instructions as title.
title(message, 'FontSize', fontSize);
bgroundROI = roipoly();
close gcf;

%% Calculate uADR and nADR
%Combine User and Isodata collected masks
somaMask = redMaskPrep & somaMaskPrep;
axonMask = redMaskPrep & cellMaskPrep &~ somaMaskPrep;

%Make max intensity projections for both red and green data
greenMax1 = max(Gdata1,[],3);
redMax2 = max(Rdata2,[],3);

%Calculate background subtracted red and green images
greenMax1 = greenMax1 - mean(greenMax1(bgroundROI));
redMax2 = redMax2 - mean(redMax2(bgroundROI));

%Calculating uADR and ADR for red probe
D_g = mean(greenMax1(somaMask));
A_g = mean(greenMax1(axonMask));
uADR = A_g/D_g;

D_r = mean(redMax2(somaMask));
A_r = mean(redMax2(axonMask));
redADR = A_r/D_r;

%Normalized ADR
nADR = uADR/redADR;

%Outputting number of pixels in axon as proxy for arbor size
axonSize = sum(axonMask(:));

%% Make concentric cirle mask over
%User selects soma center and concentric circle mask is made
circlemask = ConcentricCircleMaker(redMax1,radius);

%Get axon intensity values as a function of radial distance from soma
%center
[avg_v_radG,avg_v_radR,numpix,axConCir]...
    = AvgFromCircleMask(circlemask,axonMask,greenMax1,redMax2);

%Make imagesc-like picture of binned axon
axConCir(somaMask) = length(numpix) + 1;
%make a borderless image
iptsetpref('ImshowBorder','tight');
figure;
image(axConCir)
set(gcf,'Colormap',jet(max(axConCir(:))+1));
newCmap = colormap;
newCmap(1,:) = [0,0,0];
colormap(newCmap);
%reset border settings
iptsetpref('ImshowBorder','loose');

%Get uADR and nADR binned by distance from the soma
uADR_v_rad = [1,avg_v_radG/D_g];%1 is for circle containing soma
redADR_v_rad = [1,avg_v_radR/D_r];
nADR_v_rad = uADR_v_rad./redADR_v_rad;

%plot uADR and nADR scatters as a function of distance from soma
x_plot = 0:targetRadius:length(numpix)*targetRadius;%plot each distance from soma bin
numpix_adr = [sum(somaMask(:)),numpix];%add element for soma
%set up color map, soma color first
newCmap([3:end,2],:) = newCmap([2,3:end],:);%soma color first
newCmap(1,:) = [];%no need for black background
%plot uADR as a function of distance from soma
figure;
[axu, p1u, p2u]=plotyy(x_plot,uADR_v_rad,x_plot,numpix_adr);hold on;
scatter(x_plot,uADR_v_rad,[],newCmap,'filled');hold off;
xlabel('Distance from soma (\mum)')
set(get(axu(1),'Ylabel'),'String','Un-normalized ADR')
set(get(axu(2),'Ylabel'),'String','Pixels per bin')
%plot nADR as a function of distance from the soma
figure;
[axn, p1n, p2n]=plotyy(x_plot,nADR_v_rad,x_plot,numpix_adr);hold on;
scatter(x_plot,nADR_v_rad,[],newCmap,'filled')
xlabel('Distance from soma (\mum)')
set(get(axn(1),'Ylabel'),'String','Normalized ADR')
set(get(axn(2),'Ylabel'),'String','Pixels per bin')

%% Make final output variable
newADRout.x_plot = x_plot;
newADRout.nADR_v_rad = nADR_v_rad;
newADRout.uADR_v_rad = uADR_v_rad;
newADRout.redADR_v_rad = redADR_v_rad;
newADRout.numpix = numpix_adr;

%% Check for data duplication, if not package, if so, maybe replace
%Disabled because data does not need to be aggregated
%
% checkFile = exist('ADRdat.mat');
% %if file with ADR data has been created
% if checkFile==2
%     load('ADRdat.mat')
% else
%     ADRout = newADRout;
% end
% 
% [newFileTF,exchangeDatTF,exchangeIdx] =...
%             newStrucCheck(newADRout,ADRout);
% 
% %If file new, append new structure piece to FRAPraw and save
% if newFileTF
%     %If this is the first time data is entered to the structure, load in
%     %newFRAPraw
%     ADRout = [ADRout,newADRout];
%     
% %If old file, but user would like to exchange data
% elseif exchangeDatTF
%     ADRout(exchangeIdx) = newADRout;
% %Else do nothing with the new data
% end
% 
% %Package strucutre output (time series), output scalar values
% save('ADRdat.mat','ADRout')
% 
% out = [uADR,redADR,nADR,axonSize,laserPower488,laserPower561];

%% Make Overlay Figures
%make the image
figure
imshow(greenMax1);
showMaskAsOverlay2(0.2,axonMask,'y',[],0);
showMaskAsOverlay2(0.2,somaMask,'m',[],0);

%make a white background with the overlay opaque
figure
whiteB = ones(size(greenMax1));
imshow(whiteB);
showMaskAsOverlay2(1,axonMask,'y',[],0);
showMaskAsOverlay2(1,somaMask,'m',[],0);

