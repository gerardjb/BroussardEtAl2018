%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the bleach profile of dLGN->V1 boutons from YL 07082016 axon-G6s
%and 041518 syG6s datasets.
%Joey Broussard
%Tian Lab, UC Davis
%04/15/2018
%
%Makes two figure types:
%   1. A figure showing the fluorescence averaged across all ROIs and
%       trials with data for each ROI normalized to the first 100 frames
%       and boxcar filtered over 100 frames.
%   2. Heat maps of the first 70 ROIs from each dataset.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pull in help functions
path1 = fileparts(which('PlotNaInVivoBleach042318.m'));
addpath(genpath([path1,'/HelperFunctions']));

%Load data
load('InVivoBleaching042318.mat')
Img_rate = 2;%hertz
time_vect = 0:1/Img_rate:size(G6s_20130227_03_bleach,1)/Img_rate;
time_vect = time_vect(1:end-1);

%Smooth the data by a 10 point boxcar to make plots clearer
H = fspecial('average',[100,1]);
G6s_20130227_03_bleach_smooth = imfilter(G6s_20130227_03_bleach,...
    H,'replicate');
syG6s_419149set_bleachMat_smooth = ...
    imfilter(syG6s_419149set_bleachMat(1:length(time_vect),:),...
    H,'replicate');
axonG6s_20160708set_bleachMat_smooth =...
    imfilter(axonG6s_20160708set_bleachMat(1:length(time_vect),:),...
    H,'replicate');

%Fit to a double exponential
G6sDecay_cf = DecayFit_NonZero(mean(G6s_20130227_03_bleach_smooth,2),time_vect);
syG6sDecay_cf = DecayFit_NonZero(mean(syG6s_419149set_bleachMat_smooth,2),time_vect);
axonG6sDecay_cf = DecayFit_NonZero(mean(axonG6s_20160708set_bleachMat_smooth,2),time_vect);
%Make plottable version
G6sDecay_hat = feval(G6sDecay_cf,time_vect);
syG6sDecay_hat = feval(syG6sDecay_cf,time_vect);
axonG6sDecay_hat = feval(axonG6sDecay_cf,time_vect);

%Make the plots
figure;
%G6s plots
h_u = shadedErrorBar(time_vect,G6s_20130227_03_bleach_smooth',...
    {@mean,@(x) std(x)/sqrt(377)},...
    '-.b');
hold on;
plot(time_vect,G6sDecay_hat,'b');
%syG6s plots
h_sy = shadedErrorBar(time_vect,syG6s_419149set_bleachMat_smooth'...
    ,{@mean,@(x) std(x)/sqrt(70)},...
    '-.y');
plot(time_vect,syG6sDecay_hat,'y');
%axonG6s plots
h_axon = shadedErrorBar(time_vect,axonG6s_20160708set_bleachMat_smooth',...
    {@mean,@(x) std(x)/sqrt(354)},...
    '-.r');
plot(time_vect,axonG6sDecay_hat,'r');
%Add markings
xlabel('time (sec)')
ylabel('Normalized Fluorescence (A.U.)')

%% Doing the heatmaps for the individual ROIs
%uG6s
subG6s = G6s_20130227_03_bleach(1:2640,1:70);
avgFu = mean(subG6s,1);
[~,sortU] = sort(avgFu);
figure;
imagesc(subG6s(1:2640,sortU)',[0,3]);
title('uG6s sorted')

%syG6s
sub_syG6s = syG6s_419149set_bleachMat(1:2640,1:70);
avgFsy = mean(sub_syG6s,1);
[~,sortSy] = sort(avgFsy);
figure;
imagesc(sub_syG6s(1:2640,sortSy)',[0,3]);
title('syG6s Sorted Bleach')

%axonG6s
sub_axonG6s = axonG6s_20160708set_bleachMat(1:2640,1:70);
avgFaxon = mean(sub_axonG6s,1);
[~,sortAxon] = sort(avgFaxon);
figure;
imagesc(sub_axonG6s(1:2640,sortAxon)',[0,3]);
title('axonG6s Sorted Bleach')