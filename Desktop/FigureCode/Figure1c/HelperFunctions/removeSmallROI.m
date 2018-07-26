function [BWwO,cc_wO] = removeSmallROI(BWw,smallROISize,neighborhood,showIm)
%by: Joey Broussard
%at UC Davis on 09/08/14
%
%This function takes a BW image, removes small ROIs smaller than
%smallROISize and returns the resultant image. Connected Components
%detemined by neighborhood, which can have values of 4 or 8. showIm toggles
%making a figure showing cleaned image

%neiborhood default is 4
if nargin<3
    neighborhood = 4;
    %If showIm not provided, default is to not show
elseif nargin<4
    showIm = 0;
end
    

z = bwconncomp(BWw,neighborhood);
s = regionprops(z,'Area');
areas = [s.Area];
largeIdx = find(areas>smallROISize);%ROIs with an area under smallROIsize

%Make the zNew structure to hold ROI data above size criterion
zNew = struct;
zNew.Connectivity = neighborhood;
zNew.ImageSize = z.ImageSize;
zNew.NumObjects = length(largeIdx);

%Add all Idx lists above the size criterion
zNew.PixelIdxList = z.PixelIdxList(largeIdx);

%connected components output structure
cc_wO = zNew;
%binary mask output image
BWwO = zeros(size(BWw));
BWwO(cat(1,cc_wO.PixelIdxList{:})) = true;

%Show the cleaned up image if showIm toggled on
if showIm
    figure;
    imshow(BWwO)
end