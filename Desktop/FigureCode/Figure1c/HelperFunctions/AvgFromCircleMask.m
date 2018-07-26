function [avg_v_radG,avg_v_radR,numpix,axConCir] = AvgFromCircleMask(circlemask,axonMask,greendata,reddata)
%%This function takes a concentric circle mask (made by
%%ConcentricCircleMaker function) and the axon mask. It applies these to
%%a 3D imge to extract average image intensities withi each part of the
%%mask. It also returns number of pixels averaged in each section as a
%%control.

%Determine largest present circle by numeric value
numcircle = max(circlemask(:));

%Initialize Output variables built in for loop
avg_v_radG = zeros(1,numcircle);
avg_v_radR = zeros(1,numcircle);
numpix = zeros(1,numcircle);
axConCir = zeros(size(greendata));

for i = 1:numcircle
    %indices for the iTH annulus
    circlemaskThis = circlemask==i;
    %indices of axon contained in iTH annulus
    goodIdxs = circlemaskThis&axonMask;
    
    %build outputs
    axConCir(goodIdxs) = i;
    avg_v_radG(i) = mean(greendata(goodIdxs));
    avg_v_radR(i) = mean(reddata(goodIdxs));
    numpix(i) = sum(goodIdxs(:));
    
    i
end
%Clip off trailing circles that have no axon in them
lastNum = find(numpix,1,'last');
%dropping innermost circle because many file don't have one
avg_v_radG = avg_v_radG(2:lastNum);
avg_v_radR = avg_v_radR(2:lastNum);
numpix = numpix(2:lastNum);