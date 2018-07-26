function [I,mPix,nPix,nImage] = openTiff(fname,nChannel)
%
%This function is used to open 3D tif files where 3rd dim is time or Z. If
%two channels are present, nImage and I are output as a 2 element vector
%and structure respectively.

%If user doesn't input nChannels, default is 1
if nargin<2
    nChannel = 1;
end

%Getting Image Info
InfoImage=imfinfo(fname);
nPix=InfoImage(1).Width;
mPix=InfoImage(1).Height;
nImage=length(InfoImage);
bitdepth = InfoImage(1).BitDepth;
% if bitdepth~=8||bitdepth~=16||bitdepth~=32
%     bitdepth = 8;
% end
bitDepthStr = sprintf('uint%d',bitdepth);

%if here is based on the presence of one or two channels
if nChannel == 1
    
    %Initialize Image matrix
    I=zeros(mPix,nPix,nImage,bitDepthStr);

    %Input each image into a slice of the output matrix
    for iImage=1:nImage
       I(:,:,iImage)=imread(fname,'Index',iImage,'Info',InfoImage);
    end
    
elseif nChannel == 2
    
    %Initialize Image matrix for green and red channels
    nPerIm = nImage/2;%Note*:this assumes same number of images per channel
    I.green = zeros(mPix,nPix,nPerIm,bitDepthStr);
    I.red = zeros(mPix,nPix,nPerIm,bitDepthStr);
    
    %Input alternating images as a slice of each output matrix
    for iPerIm = 1:nPerIm
        I.green(:,:,iPerIm) = imread(fname,'Index',2*iPerIm-1,'Info',InfoImage);
        I.red(:,:,iPerIm) = imread(fname,'Index',2*iPerIm,'Info',InfoImage);
    end
    
    %Set nImage = nImage/2 for output of file;
    nImage = nPerIm;
end
