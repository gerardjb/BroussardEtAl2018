function pixSize = getPixSize(fname)
%This function uses the lsm1 structure produced by the lsminfo suite of
%tools to find the edge-wise size of pixels in microns as found in .lsm
%metadata.

[lsm1,~,~] = lsminfo(fname);
pixSize = lsm1.VoxelSizeX * 10^6;