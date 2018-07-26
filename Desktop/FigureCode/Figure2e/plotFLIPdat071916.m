function plotFLIPdat071916
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotFLIPdat071916
%Joey Broussard
%Tian Lab, UC Davis
%07/19/2016
%
%This function takes and plots the trial-averaged FLIP bleach curves.
% Returns:
%   1. The trial averaged plot of the bleach curves of the different
%   constructs
%   2. The trial averaged image aqcuisition bleaching
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%load the data structures
load('flipRaw.mat')

%Pull in help functions
path1 = fileparts(which('plotFLIPdat071916.m'));
addpath(genpath([path1,'/HelperFunctions']));

%Get variables by construct type
consts = nominal(strtok({flipRaw.id_const}));
constLev = getlevels(consts);

%colormap
cmap = varycolor(numel(constLev));

%build plot by each construct
for iConst = 1:numel(constLev)
    %pack all time series for a single construct into a matrix
    flipMat = cat(2,flipRaw(consts==constLev(iConst)).flip_norm);
    flipAvg = mean(flipMat,2);
    flipErr = std(flipMat,0,2)/sqrt(size(flipMat,2)+6);
    
    %figure for the FLIP region average time course
    figure(1);
    hold on;
    plot(t_vect(10:end)-t_vect(10),flipAvg(10:end),'Color',cmap(iConst,:))
    hold off;
    
    %figure to do the error curves plotted along with the stderr
    if iConst==1||iConst==2||iConst==10
        figure(4);
        hold on;
        shadedErrorBar(t_vect(10:end)-t_vect(10),flipAvg(10:end),flipErr(10:end),...
            {'Color',cmap(iConst,:)});
        hold off;
    end
    
    %pack all the imaged bleach curves into a matrix
    bleachMat = cat(2,flipRaw(consts==constLev(iConst)).bleachCheck_norm);
    bleachAvg = mean(bleachMat,2);
    bleachErr = std(bleachMat,0,2)/sqrt(size(bleachMat,2));
    
    %figure for the trial average bleach curve
    figure(2);
    hold on;
    plot(t_vect,bleachAvg,'Color',cmap(iConst,:));
    hold off;
    
    %calculate the FLIP with effect of imaging bleach subtracrted off
    flipMinMat = flipAvg - bleachAvg + 1;
    
    %figure for plotting FLIP with imaging bleach subtracted off
    figure(3);
    hold on;
    plot(t_vect,flipMinMat,'Color',cmap(iConst,:));
    
    figure(iConst + 4);
    plot(t_vect,flipMat);
    
end
%add onto the figures
figure(1);
legend(cellstr(constLev));

figure(2);
legend(cellstr(constLev));

figure(4);
legend(cellstr(constLev([1,2,10])));