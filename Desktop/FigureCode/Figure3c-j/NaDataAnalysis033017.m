function NaDataAnalysis033017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NaDataAnalysis033017
%Joey Broussard
%Tian Lab, UC Davis
%03/15/2017
%
%Processing data from the dLGN dataset.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%Pull in help functions
path1 = fileparts(which('NaDataAnalysis033017.m'));
addpath(genpath([path1,'/HelperFunctions']));

%Loading in the files
load('axonG6s_OS_boutons_0_100um.mat')%Gap data, 0-100um
load('cytoG6s_OS_boutons_0_100um.mat')%cyto data, 0-100um

%Prep the Gap data
clear gap_F0_allROI

%Prep the cyto data
nAS = size(sorted_ORI_DFF00_100_AS,4);
nDS = size(sorted_ORI_DFF00_100_DS,4);
cyto_dff0_allROI = cat(4,sorted_ORI_DFF00_100_AS,sorted_ORI_DFF00_100_DS);
clear sorted_ORI_DFF00_100_AS sorted_ORI_DFF00_100_DS

%Outputs
gapOSI = [];
gapDSI = [];
gapPrefDir = [];
gapPrefOri = [];
gapDifStd = [];%Difference between baseline and stimulus standard deviation
gapF0 = []; %median of the F values for each ROI data matrix
gapDFF_pref = []; %Averaged across trials for preferred direction
gapDFF_neg = []; %Avergaged across trials for most negative direction
gapSNR_pref = []; %gapDFF_pref divided by the std of the response period
gapSNR_neg = [];
gapMeanDFF_pref = []; %Averaged across trials for preferred direction, then across stim presentation time
gapSNR = []; %point-wise std across trials, then averaged across stim time
gapNegResp = [];%holds the value of the most negative response
gapNegRespSNR = [];%most negative response divided by the stdBase or similar
gapID = []; %hold animal ID and session number


ctlOSI = [];
ctlDSI = [];
ctlSNR = [];
ctlPrefDir = [];
ctlPrefOri = [];
ctlDifStd = [];
ctlF0 = [];
ctlDFF_pref = [];
ctlDFF_neg = [];
ctlSNR_pref = [];
ctlSNR_neg = [];
ctlMeanDFF_pref = [];
ctlSNR = [];
ctlNegResp = [];
ctlNegRespSNR = [];
ctlID = [];

%For the loop
plotON = 0;



for iConst = 1:2
    
    %Deterimne the number of ROIs in this session
    if iConst == 1
       nROI = size(gap_dff0_allROI,4);
       nDirs = size(gap_dff0_allROI,2);
    elseif iConst == 2
       nROI = size(cyto_dff0_allROI,4);
       nDirs = size(cyto_dff0_allROI,2);
    end
    
       for iROI = 1:nROI
           iROI
            
           %make vector of directions
            dirs = 0:360/nDirs:360;
            dirs = dirs(1:end-1);
           
            if iConst==1
                %Getting the dff and with the correct session and ROI number
                thisDFF = gap_dff0_allROI(:,:,:,iROI);
                
            else
                %Getting the dff with the correct session and ROI number
                thisDFF = cyto_dff0_allROI(:,:,:,iROI);

            end

            %% Calculating the oritientation and direction values and preferences
            %Function for separating DFF matrix into direction-specific
            %trials matrix and extracting average responses to each
            %dirDFF is of form trials x time x stimDir
            [dirDFF,dirResp,dirs] = calcDirResp(thisDFF,dirs);
            
            %Make a plot of the trial average responses to each stim direction for the
            %selected ROI. Toggle on and off with plotON
            if plotON
                checkRespPlot(dirs,dirResp,dirDFF,thisDFF)
            end%plotOn
            
            %Correcting for differences in acquisition settings
            gain = 1;
            
            %Calculate the indices for direction and orientation
            [prefResp,prefDirIdx] = max(dirResp);
            [~,lowDirIdx] = min(dirResp); %direction of smallest response
            [~,bestDirIdx] = max(abs(dirResp));
            [negResp,negDirIdx] = min(dirResp);
            %This is where to build the new approach to preferred
            %orientation calculation
%             if dirResp(prefDirIdx)<0
%                 dirResp = -dirResp;
%                 dirResp = dirResp - min(dirResp);
%           %Getting the step size for finding ortho and null directions
            quarterTurn = (length(dirs)-1)/4;
            prefDirIdx = prefDirIdx(1);
            orthoIdx1 = mod(prefDirIdx + quarterTurn - 1, length(dirs)-1) + 1;
            orthoIdx2 = mod(prefDirIdx - quarterTurn - 1, length(dirs)-1) + 1;
            orthoResp = mean([dirResp(orthoIdx1),dirResp(orthoIdx2)]);
            %orthoResp = dirResp(orthoIdx1);
            nullIdx = mod(prefDirIdx + 2*quarterTurn - 1, length(dirs)-1) + 1;
            nullResp = dirResp(nullIdx);
            %I turned off the zeroing out of the values because some of the
            %cells clearly do have inhibition at some of the orientations
%             if nullResp<0
%                 nullResp = 0;
%             elseif orthoResp<0
%                 orthoResp = 0;
%             end
            %Calculate OSI and DSI
            thisOSI = (prefResp - orthoResp)/(prefResp + orthoResp);
            thisDSI = (prefResp - nullResp)/(prefResp + nullResp);
            thisPrefDir = dirs(prefDirIdx);
            thisPrefOri = dirs(mod(prefDirIdx-1,(length(dirs)-1)/2) + 1);
                       
            %% Calculating F0, noise, and SNR parameters
            %Calculate F, DFF_prefered, and SNR
            %DFF_preferred is the response in the stim period
            thisDFF_pref = squeeze(nanmean(dirDFF(:,:,prefDirIdx),1));
            thisMeanDFF_pref = mean(thisDFF_pref(4:11));
            thisDFF_neg = squeeze(nanmean(dirDFF(:,:,negDirIdx),1));
            %SNR calculations
            std_vect = squeeze(nanstd(dirDFF(:,:,bestDirIdx),1));
            std_vect_low = squeeze(nanstd(dirDFF(:,:,lowDirIdx(1)),1));
            std_vect_neg = squeeze(nanstd(dirDFF(:,:,negDirIdx),1));
            stdBase = mean(std_vect(1:3));%4:11));
            stdBaseNeg = mean(std_vect(1:3));%4:11));
            std_denom = mean(std_vect_low);%4:11));%per Leo's request
            std_dif = std_denom - stdBase;
            thisSNR_pref = thisDFF_pref/std_denom;
            thisSNR_neg = thisDFF_neg/stdBaseNeg;
            thisSNR = thisMeanDFF_pref/std_denom;
            thisNegRespSNR = negResp/std_denom;
            
%             %% Statistical significance tests
%             thisDFFbase = thisDFF_pref(1:round(3));
%             thisDFFresp = thisDFF_pref(5:11);
%             thisID(4) = ttest(thisDFFbase,thisDFFresp,'Alpha',0.01);
            
            
           %% Packaging output variables
           %Packaging the ROI-specific outputs
           if iConst==1
                gapOSI = cat(1,gapOSI,thisOSI);
                gapDSI = cat(1,gapDSI,thisDSI);
                gapPrefDir = cat(1,gapPrefDir,thisPrefDir);
                gapPrefOri = cat(1,gapPrefOri,thisPrefOri);
                gapDifStd = cat(1,gapDifStd,std_dif);
                gapDFF_pref = nancat(1,gapDFF_pref,thisDFF_pref); %Averaged across trials, then across stim time
                gapSNR_pref = nancat(1,gapSNR_pref,thisSNR_pref);
                gapDFF_neg = nancat(1,gapDFF_neg,thisDFF_neg);
                gapSNR_neg = nancat(1,gapSNR_neg,thisSNR_neg);
                gapMeanDFF_pref = cat(1,gapMeanDFF_pref,thisMeanDFF_pref); 
                gapSNR = cat(1,gapSNR,thisSNR); %point-wise std across trials, then averaged across stim time
                gapNegResp = cat(1,gapNegResp,negResp);
                gapNegRespSNR = cat(1,gapNegRespSNR,thisNegRespSNR);
           elseif iConst==2
                ctlOSI = cat(1,ctlOSI,thisOSI);
                ctlDSI = cat(1,ctlDSI,thisDSI);
                ctlPrefDir = cat(1,ctlPrefDir,thisPrefDir);
                ctlPrefOri = cat(1,ctlPrefOri,thisPrefOri);
                ctlDifStd = cat(1,ctlDifStd,std_dif);
                ctlDFF_pref = nancat(1,ctlDFF_pref,thisDFF_pref); %Averaged across trials, then across stim time
                ctlDFF_neg = nancat(1,ctlDFF_neg,thisDFF_neg);
                ctlSNR_pref = nancat(1,ctlSNR_pref,thisSNR_pref);
                ctlSNR_neg = nancat(1,ctlSNR_neg,thisSNR_neg);
                ctlMeanDFF_pref = cat(1,ctlMeanDFF_pref,thisMeanDFF_pref);
                ctlSNR = cat(1,ctlSNR,thisSNR); %point-wise std across trials, then averaged across stim time                                
                ctlNegResp = cat(1,ctlNegResp,negResp);
                ctlNegRespSNR = cat(1,ctlNegRespSNR,thisNegRespSNR);
           end
           
       end%iROI
       
end%iConst

%% output graphs as ECDFs
%Load and normalize the F0 data
load('unormalizedF0values.mat')
gapNormF0 = gap_G6s_F0(:,5:end)./repmat(gap_G6s_F0(:,4),1,size(gap_G6s_F0,2)-4).^2;
ctlNormF0 = ctl_G6s_F0(:,5:end)./repmat(ctl_G6s_F0(:,4),1,size(ctl_G6s_F0,2)-4).^2;

%ECDF for F0, DFF_pref, SNR
figure;
subplot(2,2,1)
plotCDF(ctlNormF0(:),'b');hold on;
plot([median(ctlNormF0(:)),median(ctlNormF0(:))],[0,1],'b');%add on median line
text(median(ctlNormF0(:)),.25,num2str(median(ctlNormF0(:))),'color',[0,0,1]);%and text
plotCDF(gapNormF0(:),'r');
plot([median(gapNormF0(:)),median(gapNormF0(:))],[0,1],'r');%add on median line
text(median(gapNormF0(:)),.75,num2str(median(gapNormF0(:))),'color',[1,0,0])%and text
title('F0')
subplot(2,2,2)
plotCDF(ctlMeanDFF_pref,'b');
hold on;
plot([median(ctlMeanDFF_pref),median(ctlMeanDFF_pref)],[0,1],'b');%add on median line
text(median(ctlMeanDFF_pref),.25,num2str(median(ctlMeanDFF_pref)),'color',[0,0,1]);%and text
plotCDF(gapMeanDFF_pref,'r');
plot([median(gapMeanDFF_pref),median(gapMeanDFF_pref)],[0,1],'r');%add on median line
text(median(gapMeanDFF_pref),.75,num2str(median(gapMeanDFF_pref)),'color',[1,0,0])%and text
title('Preferred DFF')
subplot(2,2,3)
plotCDF(ctlSNR,'b');hold on;
plot([median(ctlSNR),median(ctlSNR)],[0,1],'b');%add on median line
text(median(ctlSNR),.25,num2str(median(ctlSNR)),'color',[0,0,1]);%and text
plotCDF(gapSNR,'r');
plot([median(gain*gapSNR),median(gain*gapSNR)],[0,1],'r');%add on median line
text(median(gain*gapSNR),.75,num2str(median(gain*gapSNR)),'color',[1,0,0])%and text
title('SNR')
legend('ctl','gap')
subplot(2,2,4)
% plotCDF(ctl

% %CDF for OSI and DSI
% figure;
% subplot(1,2,1)
% plotCDF(ctlOSI,'b');hold on;
% plotCDF(gapOSI,'r');
% title('OSI')
% subplot(1,2,2)
% plotCDF(ctlDSI,'b');hold on;
% plotCDF(gapDSI,'r');
% title('DSI')
% legend('ctl','gap')

%% Individual ROI traces
%Individual ROI preferred direction DFF/SNR trace maps
%DFF preferred
makeDFF_SNRheatmap(gapDFF_pref,ctlDFF_pref,'preferred DFF')
%SNR preferred
makeDFF_SNRheatmap(gain*gapSNR_pref,ctlSNR_pref,'preferred SNR')
%DFF most negativex`
makeDFF_SNRheatmap(gapDFF_neg,ctlDFF_neg,'most negative DFF')
%SNR most negative
makeDFF_SNRheatmap(gapSNR_neg,ctlSNR_neg,'most negative SNR')

%% Preferred direction and orientation graphs
%preferred direction
figure;
histogram(ctlPrefDir,-22.5:45:337.5,'Normalization','PDF');
hold on;histogram(gapPrefDir,-22.5:45:337.5,'Normalization','PDF');
set(gca,'xtick',0:45:315);
legend('ctl','gap')
title('Preferred Direction')

%preferred orientation
figure;
histogram(ctlPrefOri,-22.5:45:157.5,'Normalization','PDF');
hold on;histogram(gapPrefOri,-22.5:45:337.5,'Normalization','PDF');
set(gca,'xtick',0:45:135);
legend('ctl','gap')
title('Preferred Orientation')

% %Broken out by animals
% %Gap preferred direction
% breakByAnimal(gapPrefDir,-22.5:45:337.5,gapID(:,2),'gap');
% %Ctl preferred direction
% breakByAnimal(ctlPrefDir,-22.5:45:337.5,ctlID(:,2),'ctl');
% %Gap preferred orientation
% breakByAnimal(gapPrefOri,-22.5:45:157.5,gapID(:,2),'gap');
% %Ctl preferred orientation
% breakByAnimal(ctlPrefOri,-22.5:45:157.5,ctlID(:,2),'ctl');

% %% Breaking out data by animal/session
% %DSI by session
% figure;
% for i = 1:nSession(1)
%     if i==1;
%         [~,~,h1] = plotCDF(gapDSI(gapID(:,1)==i),'r');
%     else
%         [~,~,~] = plotCDF(gapDSI(gapID(:,1)==i),'r');
%     end
%     hold on;
% end
% hold on;
% for i = 1:nSession(2)
%     if i==1
%         [~,~,h2] = plotCDF(ctlDSI(ctlID(:,1)==i),'b');
%     else
%         plotCDF(ctlDSI(ctlID(:,1)==i),'b');
%     end
%     hold on
% end
% legend([h1,h2],{'gap','ctl'})
% xlim([0,1])
% ylim([0,1])
% xlabel('DSI')
% ylabel('Cumulative probability')
% title('DSI broken out by session')
% 
% %OSI broken out by session
% figure;
% for i = 1:nSession(1)
%     if i==1;
%         [~,~,h1] = plotCDF(gapOSI(gapID(:,1)==i),'r');
%     else
%         [~,~,~] = plotCDF(gapOSI(gapID(:,1)==i),'r');
%     end
%     hold on;
% end
% hold on;
% for i = 1:nSession(2)
%     if i==1
%         [~,~,h2] = plotCDF(ctlOSI(ctlID(:,1)==i),'b');
%     else
%         plotCDF(ctlOSI(ctlID(:,1)==i),'b');
%     end
%     hold on
% end
% legend([h1,h2],{'gap','ctl'})
% xlim([0,1])
% ylim([0,1])
% xlabel('OSI')
% ylabel('Cumulative probability')
% title('OSI broken out by session')
% 
% %DSI by animal
% figure;
% for i = 1:max(gapAnimalID)
%     if i==1;
%         [~,~,h1] = plotCDF(gapDSI(gapID(:,2)==i),'r');
%     else
%         [~,~,~] = plotCDF(gapDSI(gapID(:,2)==i),'r');
%     end
%     hold on;
% end
% hold on;
% for i = 1:max(ctlAnimalID)
%     if i==1
%         [~,~,h2] = plotCDF(ctlDSI(ctlID(:,2)==i),'b');
%     else
%         plotCDF(ctlDSI(ctlID(:,2)==i),'b');
%     end
%     hold on
% end
% legend([h1,h2],{'gap','ctl'})
% xlim([0,1])
% ylim([0,1])
% xlabel('DSI')
% ylabel('Cumulative probability')
% title('DSI broken out by animal')
% 
% %OSI broken out by animal
% figure;
% for i = 1:max(gapAnimalID)
%     if i==1;
%         [~,~,h1] = plotCDF(gapOSI(gapID(:,1)==i),'r');
%     else
%         [~,~,~] = plotCDF(gapOSI(gapID(:,1)==i),'r');
%     end
%     hold on;
% end
% hold on;
% for i = 1:max(ctlAnimalID)
%     if i==1
%         [~,~,h2] = plotCDF(ctlOSI(ctlID(:,2)==i),'b');
%     else
%         plotCDF(ctlOSI(ctlID(:,2)==i),'b');
%     end
%     hold on
% end
% legend([h1,h2],{'gap','ctl'})
% xlim([0,1])
% ylim([0,1])
% xlabel('OSI')
% ylabel('Cumulative probability')
% title('OSI broken out by animal')

% %% Binning OSI and DSI by SNR
% %OSI with gap binned
% binOSI_DSIbySNR(gapOSI,gapSNR,ctlOSI,'OSI','gap')
% %OSI with ctl binned
% binOSI_DSIbySNR(ctlOSI,ctlSNR,gapOSI,'OSI','ctl')
% %DSI with gap binned
% binOSI_DSIbySNR(gapDSI,gapSNR,ctlDSI,'DSI','gap')
% %DSI with ctl binned
% binOSI_DSIbySNR(ctlDSI,ctlSNR,gapDSI,'DSI','ctl')

% %% SNR binning direct comparison across contructs
% nBin = 5;
% 
% figure;
% cmap1 = winter(nBin);
% cmap2 = spring(nBin);
% legStr1 = [];
% legStr2 = [];
% legHand1 = [];
% legHand2 = [];
% for i = 1:nBin
%     redmap = [1 - 0.5*(1-i/nBin),0.4*(1-i/nBin),0.4*(1-i/nBin)];
%     tholdLow = quantile(gapSNR,1/nBin*i - 1/nBin);
%     tHoldHigh = quantile(gapSNR,1/nBin*i);
%     [~,~,h1] = plotCDF(gapOSI(gapSNR>tholdLow & gapSNR<=tHoldHigh),redmap);
%     legStr1 = cat(1,legStr1,{[num2str(tholdLow),' ',num2str(tHoldHigh),' gap']});
%     legHand1 = cat(1,legHand1,h1);
%     hold on;
%     
%     bluemap = [0.4*(1-i/nBin),0.4*(1-i/nBin),1 - 0.5*(1-i/nBin)];
%     tholdLow = quantile(ctlSNR,1/nBin*i - 1/nBin);
%     tHoldHigh = quantile(ctlSNR,1/nBin*i);
%     [~,~,h2] = plotCDF(ctlOSI(ctlSNR>tholdLow & ctlSNR<=tHoldHigh),bluemap);
%     legStr2 = cat(1,legStr2,{[num2str(tholdLow),' ',num2str(tHoldHigh),' ctl']});
%     legHand2 = cat(1,legHand2,h2);
% end
% 
% xlim([0,1])
% xlabel('OSI')
% ylabel('Cumulative Probabilty')
% legend([legHand1;legHand2],[legStr1;legStr2]);
% title('Binning ROIs by SNR values')

% %% Scatters of SNR agaisnt OSI, scatter hist of Negresp and SNR
% %Scatterhist for negative response and SNR, DFF based
% figure;
% scatterhist([ctlNegResp(ctlNegResp<0);gapNegResp(gapNegResp<0)],...%negative responses below 0
%     [ctlSNR(ctlNegResp<0);gapSNR(gapNegResp<0)],...%corresponding SNRs
%     'Group',[zeros(length(ctlSNR(ctlNegResp<0)),1);ones(length(gapSNR(gapNegResp<0)),1)]);%grouping
% xlabel('Negative Response')
% ylabel('SNR')
% legend('gap','ctl')
% title('Scatterhist of SNR and negative responses')
% 
% %Scatterhist for negative response and SNR, SNR based
% figure;
% scatterhist([ctlNegRespSNR(ctlNegRespSNR<0);gapNegRespSNR(gapNegRespSNR<0)],...%negative responses below 0
%     [ctlSNR(ctlNegRespSNR<0);gapSNR(gapNegRespSNR<0)],...%corresponding SNRs
%     'Group',[zeros(length(ctlSNR(ctlNegResp<0)),1);ones(length(gapSNR(gapNegResp<0)),1)]);%grouping
% xlabel('Negative Response')
% ylabel('SNR')
% legend('gap','ctl')
% title('Scatterhist of SNR and negative responses')

%% Function Block

function [dirDFF,dirResp,dirs] = calcDirResp(thisDFF,dirs)
%Separates DFF matrix into a 3D matrix of form trials x time x stimDir
%Separating the responses based on the direction
dirDFF = permute(thisDFF,[3,1,2]);

%Break out the maximum response from each direction
dirMean = squeeze(nanmean(dirDFF(:,:,:),1));%Note: tried averageing across the time domain first, then across trials. But to fit with the std calculations, I did it across the trials then time insteasd
dirResp = mean(dirMean(4:end,:),1);%Note: per Na, stim applied from frame 4 to 11

%Append first value to end to ease modulo calculations and
%plots
dirResp(length(dirResp) + 1) = dirResp(1);%Add on the 0 response for 360 to complete circle
dirs(length(dirs) + 1) = 360;%add on the 360 to complete the circle

function checkRespPlot(dirs,dirResp,dirDFF,thisDFF)
%Polar plot of the response as a function of stimulus orientation
figure;
polar(radians(dirs),dirResp);

%Plot the data all lined up
figure;
plot(reshape(thisDFF,size(thisDFF,1)*size(thisDFF,2)*size(thisDFF,3),1))

%avgPlots of the responses at each orientation
f1 = figure;
for i = 1:length(dirs)-1
    %input this avgplot to the figure
    subplot(1,length(dirs-1),i);
    avgPlot(dirDFF(:,:,i)');
    ylim([min(thisDFF(:)),max(thisDFF(:))]);
    thisDir = num2str(dirs(i));
    xlabel([thisDir,'^o'])
    if i==1
        ylabel('\DeltaF/F')
    end
end

%avgPlot of the trial averages with binning by stimulus orientation
%Effectively, this just gives a global sense of how reponsive this bouton
%wwas
figure;
avgPlot(squeeze(nanmean(dirDFF(:,:,:),1)))

function makeDFF_SNRheatmap(gap,ctl,descriptor)
%This function makes a heatmap of the input traces for gap and ctl. The
%heatmap varies from red to green with 0 set as black. descriptor describes
%what data is displayed, e.g. 'preferred DFF', 'most negative SNR'

%Individual ROI preferred direction DFF trace maps
%Process the data
%Convert NaN values to zeros
gap(isnan(gap)) = 0;
ctl(isnan(ctl)) = 0;
%align by maximum deflection of DFF
maxGap = mean(gap(:,4:11),2);
[~,gapSortIdx] = sort(maxGap,'Descend');
maxCtl = mean(ctl(:,4:11),2);
[~,ctlSortIdx] = sort(maxCtl,'Descend');
%color limits based on the higher of the mins and maxes from both matrices
clims = [min(ctl(:)),max(ctl(:))];
clims = [min([min(ctl(:)),min(gap(:))]),max([max(ctl(:)),max(gap(:))])];
bluePercent = abs(clims(1))/(clims(2) + abs(clims(1)));
redPercent = clims(2)/(clims(2) + abs(clims(1)));
%Stretch the colorspace over the clims with zero black
redCol = [linspace(0,1,floor(bluePercent*256)),ones(1,ceil(256*redPercent))];
greenCol = [linspace(0,1,floor(bluePercent*256)),linspace(1,0,ceil(256*redPercent))];
blueCol = [ones(1,floor(bluePercent*256)),linspace(1,0,ceil(256*redPercent))];
RWBcolormap = [redCol;greenCol;blueCol]';
figure;
%Gap subplot
subplot(1,2,1)
imagesc(gap(gapSortIdx,:));
caxis(clims);
colormap(RWBcolormap);
colorbar;
xlabel('time (A.U.)');
ylabel('ROI number');
title('GAP43-GCaMP6s')
%ctl subplot
subplot(1,2,2)
imagesc(ctl(ctlSortIdx,:));
caxis(clims);
colormap(RWBcolormap);
colorbar;
xlabel('time (A.U.)');
ylabel('ROI number');
title('GCaMP6s')
suptitle(['Preferred direction selected by maximum response:',descriptor])

function binOSI_DSIbySNR(binnedGroup,binningSNR,comparisonGroup,O_or_D,binnedConst)
%This function plots CDFs of OSI or DSI from the binnedGroup broken out by
%their quantile (currently done by 10% binning) as compared to the
%comparison group CDF. CDFs for the binned groups are on Parula, the
%comparison group is in red. O_or_D and binnedConst are char descriptors
%for OSI or DSI and the construct being binned, respectively.

%Parse which construct is being binned, the other construct is then the
%comparison
if strcmp (binnedConst,'gap')
    compConst = 'ctl';
else
    compConst = 'gap';
end


figure;
cmap = parula(10);
for i = 1:10
    tholdLow = quantile(binningSNR,0.1*i - 0.1);
    tHoldHigh = quantile(binningSNR,0.1*i);
    plotCDF(binnedGroup(binningSNR>tholdLow & binningSNR<=tHoldHigh),cmap(i,:));
    hold on;
end
%Append ctl OSI
plotCDF(comparisonGroup,'r');
xlim([0,1])
xlabel(O_or_D)
ylabel('Cumulative Probabilty')
legend('0-10% SNR','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%',...
    '70-80%','80-90%','90-100%',[compConst,O_or_D]);
title(['Binning ',O_or_D,' ROIs from ',binnedConst,' dataset'])
