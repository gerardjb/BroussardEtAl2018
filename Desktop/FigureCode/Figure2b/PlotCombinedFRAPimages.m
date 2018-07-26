%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make combined FRAP averaged line image
%Joey Broussard
%Tian Lab, UC Davis
%07/12/2017
%
%This script makes an imagesc of the FRAP data from the G6m, GAP-G6m, and
%syph-G6m experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%loading in the data
load G6mcell3_dffImg.mat
G6mDFF_1avg = dff_1avg;
load GAPcell3_dffImg.mat
GapDFF_1avg = dff_1avg;
load syph-G6cell3_dffImg.mat
syphDFF_1avg = dff_1avg;

%Making the colormap limits for the images
allImg = cat(3,G6mDFF_1avg,GapDFF_1avg,syphDFF_1avg)+1;
cmap = [min(allImg(:)),max(allImg(:))];

%reading in the pixel size and time bin length
load scale.mat

figure;
%G6m
subplot(1,3,1)
imagesc(G6mDFF_1avg+1,cmap);
%set axes
set(gca,'Yticklabel',{[50:50:400]*tScale})
set(gca,'Xticklabel',{[5:5:45]*xScale})
colorbar
%Gap-G6m
subplot(1,3,2)
imagesc(GapDFF_1avg+1,cmap);
%set axes
set(gca,'Yticklabel',{[50:50:400]*tScale})
set(gca,'Xticklabel',{[5:5:45]*xScale})
colorbar
%Gap-G6m
subplot(1,3,3)
imagesc(syphDFF_1avg+1,cmap);
%set axes
set(gca,'Yticklabel',{[50:50:400]*tScale})
set(gca,'Xticklabel',{[5:5:45]*xScale})
colorbar