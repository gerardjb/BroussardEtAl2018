%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Making the plots of Fluorescence normalized to laser power as a function
%of depth for the targeted and untargeted GCaMPs
%Joey Broussard
%Tian Lab, UC Davis
%03/15/2017
%
%Structure of variables:
%   col 1: date identifier
%   col 2: binned depth, 0<=d<100, 100<=d<200, etc.
%   col 3: depth per file name
%   col 4: laser power per file name
%   col 5-end: average F from ROIs from brightest image features
%
%Makes plots and analyses:
%   1 - figures with normalized fluorescence as a function of depth
%   2 - ANCOVA to compare the difference in fluorescence profiles as a
%       function of depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%open data
load('normFvalues080817.mat')

%Plotting normalized F as a function of depth
nConst = 2;
for iConst = 1:nConst
    %Select the martrix for this construct
    if iConst==1
        thisMat = G6s_values;
    else
        thisMat = GAP_values;
    end
    
    %Get unique depth vector for averaging purposes
    uniqueDepth = unique(thisMat(:,2));
    nDepth = length(uniqueDepth);
    nPoints = size(thisMat(:,5:end),2);
    
    %Loop over the depth values to average and make a new vector with error
    %bars, etc.
    %Receiving variables
    for iDepth = 1:nDepth
        %get correct rows containing data for each depth
        goodIdx = thisMat(:,2)==uniqueDepth(iDepth);
        
        %Prep data for averaging/error
        thisNormF_prep = thisMat(goodIdx,5:end)./repmat(thisMat(goodIdx,4),1,nPoints).^2;
        
        %Get mean and SEM calculating denominator as number of sessions
        thisNormF(iDepth) = mean(thisNormF_prep(:));
        thisNormSEM(iDepth) = std(thisNormF_prep(:))/sqrt(sum(goodIdx));
    end%iDepth
    
    %Output to appropriate variable
    if iConst==1
        G6s_normF = thisNormF;
        G6s_SEM = thisNormSEM;
        G6s_depth = uniqueDepth;
    else
        GAP_normF = thisNormF;
        GAP_SEM = thisNormSEM;
        GAP_depth = uniqueDepth;
    end%iConst
    
end%iConst

%Doing a log-linear regression to the F0 data
%Make and re-shape the data to vectors
G6s_depth_all = repmat(G6s_values(:,2),1,nPoints);
G6s_normF_all = G6s_values(:,5:end)./repmat(G6s_values(:,4),1,nPoints).^2;
G6s_x = [ones(length(G6s_depth_all(:)),1),G6s_depth_all(:)];
G6s_y = log10(G6s_normF_all(:));
GAP_depth_all = repmat(GAP_values(:,2),1,nPoints);
GAP_normF_all = GAP_values(:,5:end)./repmat(GAP_values(:,4),1,nPoints).^2;
GAP_x = [ones(length(GAP_depth_all(:)),1),GAP_depth_all(:)];
GAP_y = log10(GAP_normF_all(:));
%Do regressions
c_G6s = G6s_x\G6s_y;
c_GAP = GAP_x\GAP_y;
yHat_G6s = G6s_x*c_G6s;
yHat_GAP = GAP_x*c_GAP;

%Make the plots
%G6s
figure;
%datapoints
h1 = plot(G6s_values(:,2),...
    G6s_values(:,5:end)./repmat(G6s_values(:,4),1,nPoints).^2,'bo');
set(gca,'Yscale','log')
%Population summaries
hold on;h2 = errorbar(G6s_depth,G6s_normF,G6s_SEM,...
    'LineStyle','None',...
    'Marker','o',...
    'Color',[0,0,0.5]);
%Linear regressions
axis([0 500 0.01 100])
xlims = xlim;
h3 = plot(xlims,10.^([1,xlims(1);1,xlims(2)]*c_G6s),'b');
%legend and pretty
legend([h1(1),h2,h3],'Binned values','Bin means','Log-linear fit')
box off
grid on

%GAP
figure;
%datapoints
h4 = plot(GAP_values(:,2),...
    GAP_values(:,5:end)./repmat(GAP_values(:,4),1,nPoints).^2,'ro');
set(gca,'Yscale','log')
%Population summaries
hold on;h5 = errorbar(GAP_depth,GAP_normF,GAP_SEM,...
    'LineStyle','None',...
    'Marker','o',...
    'Color',[0.5,0,0]);
%Linear regressions
xlims = xlim;
h6 = plot(xlims,10.^([1,xlims(1);1,xlims(2)]*c_GAP),'r');
axis([0 500 0.01 100])
%legend and pretty
legend([h4(1),h5,h6],'Binned values','Bin means','Log-linear fit')
box off
grid on

%Stats - Ancova on the log-transformed data
%Setting up the data for use in the aoctool
x = [G6s_x(:,2);GAP_x(:,2)];
y = [G6s_y;GAP_y];
g = [repmat({'G6s'},[length(G6s_y),1]);repmat({'GAP'},[length(GAP_y),1])];
%Run stats
aoctool(x,y,g)