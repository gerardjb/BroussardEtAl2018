function [laserP] = getLaser(fname)
%This function gets the laser power out of the metadata unpacked by the lsm
%suite of tools. It takes the first output of the suite (a structure with a
%large amount of metadata) and returns the laser power.

[lsm1,~,~] = lsminfo(fname);
laserP = lsm1.ScanInfo.POWER{1};