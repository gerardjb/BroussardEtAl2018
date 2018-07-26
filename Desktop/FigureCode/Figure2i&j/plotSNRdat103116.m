function plotSNRdat103116
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotSNRdat103116
%Joey Broussard
%Tian Lab, UC Davis
%10/31/2016
%
%This function takes and plots data as packaged into structures with
%organization as those produced by processSNR.m. Currently returns the
%following plots:
%   1. DFF for the all constructs on each graph. One graph per FP number
%   2. Laser power bargraph.
%   3. DFF as a function of FP number per construct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%Load structure with data
load('SNRraw1.mat');load('SNRraw2.mat');load('SNRraw3.mat');load('SNRraw4.mat');
SNRraw = [SNRraw1,SNRraw2,SNRraw3,SNRraw4];clear SNRraw1 SNRraw2 SNRraw3 SNRraw4

%Pull in help functions
path1 = fileparts(which('plotSNRdat103116.m'));
addpath(genpath([path1,'/HelperFunctions']));

%Get variables for construct, stim number
%% Construct
%maybe reorder the levels to make for a better graphs
consts = nominal(strtok({SNRraw.id_const}));
constLev = getlevels(consts);
%Stimulus
stimOrd = ordinal(...
    cellfun(@num2str,{SNRraw(1).trial.id_FPnum}...
    ,'UniformOutput',false));
stimOrd = reorderlevels(stimOrd,{'1','2','5','10','20','40','80'});
stimLev = getlevels(stimOrd);
%RoiType
roiType = nominal(strtok({SNRraw.id_ROI}));
roiLev = getlevels(roiType);
%% 

%Get times to look for maximum DFF in each trace
time_vect = SNRraw(1).trial(1).time_vect;
goodTimesTF = time_vect>5&time_vect<7.5;

%Set up color scheme for use in DFF plots, initialize variable for
%RedF and ROIsizes
colConst = lines(numel(constLev));
RedF_all = [];
RedF_split = [];
ROIsizes_all = [];
ROIsizes_split = [];
t_level = 10;%threshold for laser power to image

%Prep DFF data, do plots for averaged time series
for iStim = 1:numel(stimLev)
    
    %Iitialize matrix to hold legend handle
    legHandDFF_ts = [];
    legHandSNR_ts = [];
    legHandDFF_sc = [];
    legHandSNR_sc = [];
    legHandHist = [];
    
    %Initialize matrix to hold the DFF for the histogram
    DFF_maxHist_split = [];
    SNR_maxHist_split = [];
    F0_mat_split = [];
    
    for iConst = 1:numel(constLev)
        
        %Get the laser power and ROIsizes for each construct into a matrix
        if iStim == 1%only take laser values first time through
            %RedF
            RedF_all = nancat(2,...
                RedF_all,...
                [SNRraw(roiType==roiLev(1)&consts==constLev(iConst))...
                .RedF]');
            RedF_split = nancat(3,...
                RedF_split,...
                nancat(2,...
                SNRraw(roiType==roiLev(2)&consts==constLev(iConst))...
                .RedF));
            %Areas
            ROIsizes_all = nancat(2,...
                ROIsizes_all,...
                [SNRraw(roiType==roiLev(1)&consts==constLev(iConst))...
                .Area]');
            ROIsizes_split = nancat(3,...
                ROIsizes_split,...
                nancat(2,...
                SNRraw(roiType==roiLev(2)&consts==constLev(iConst))...
                .Area));
            %Establish number of ROIs
            nROIs(iConst) = sum(~isnan(RedF_split(:))) + 1;
        end%iStim == 1
        
        
        %Prep the structure for extracting relevant data
        SNRprep_all = [SNRraw(...
            roiType==roiLev(1)&...
            consts==constLev(iConst))...
            .trial];
        stimOrd_all = ordinal(...
            cellfun(@num2str,{SNRprep_all.id_FPnum}...
            ,'UniformOutput',false));
        stimOrd_all = reorderlevels(stimOrd_all,{'1','2','5','10','20','40','80'});
        stimLev_all = getlevels(stimOrd_all);
        SNRprep_split = [SNRraw(...
            roiType==roiLev(2)&...
            consts==constLev(iConst))...
            .trial];
        stimOrd_split = ordinal(...
            cellfun(@num2str,{SNRprep_all.id_FPnum}...
            ,'UniformOutput',false));
        stimOrd_split = reorderlevels(stimOrd_split,{'1','2','5','10','20','40','80'});
        stimLev_split = getlevels(stimOrd_split);
        
        %Get the DFF and error values for each stim
        %Load in DFF, truncate off those ROIs with F0 below 2 or above 150
        F_mat_all = [SNRprep_all(stimOrd_all==stimLev_all(iStim))...
            .F];
        F_mat_split = [SNRprep_split(stimOrd_split==stimLev_split(iStim))...
            .F];
        F0_mat_all = mean(F_mat_all(1:12,:),1);
        F0_mat_split = nancat(1,F0_mat_split,mean(F_mat_split(1:12,:),1));

        %DFF
        DFF_mat_all = [SNRprep_all(stimOrd_all==stimLev_all(iStim))...
            .DFF];
        DFF_mat_split = [SNRprep_split(stimOrd_split==stimLev_split(iStim))...
            .DFF];

        %SNR
        SNR_mat_all = calcSNR(F_mat_all,1,0,12);
        SNR_mat_split = calcSNR(F_mat_split,1,0,12);
        
        %Get avg and err
        %DFF
        DFF_avg_all = mean(DFF_mat_all,2);
        DFF_err_all = std(DFF_mat_all,0,2)/sqrt(size(DFF_mat_all,2));
        DFF_avg_split = mean(DFF_mat_split,2);
        DFF_err_split = std(DFF_mat_split,0,2);%/sqrt(size(DFF_mat_split,2));
        %SNR
        SNR_avg_all = mean(SNR_mat_all,2);
        SNR_err_all = std(SNR_mat_all,0,2)/sqrt(size(DFF_mat_all,2));
        SNR_avg_split = mean(SNR_mat_split,2);
        SNR_err_split = std(DFF_mat_split,0,2);%/sqrt(size(DFF_mat_split,2));
        
        %put the values for each max DFF and the corresponding error into a
        %matrix
        %DFF
        DFF_max_all(iStim,iConst) = max(DFF_avg_all(goodTimesTF));
        DFF_errMat_all(iStim,iConst) = ...
            DFF_err_all(DFF_avg_all==DFF_max_all(iStim,iConst));
        DFF_max_split(iStim,iConst) = max(DFF_avg_split(goodTimesTF));
        DFF_errMat_split(iStim,iConst) = ...
            DFF_err_split(DFF_avg_split==DFF_max_split(iStim,iConst));
        %SNR
        SNR_max_all(iStim,iConst) = max(SNR_avg_all(goodTimesTF));
        SNR_errMat_all(iStim,iConst) = ...
            SNR_err_all(SNR_avg_all==SNR_max_all(iStim,iConst));
        SNR_max_split(iStim,iConst) = max(SNR_avg_split(goodTimesTF));
        SNR_errMat_split(iStim,iConst) = ...
            SNR_err_split(SNR_avg_split==SNR_max_split(iStim,iConst));
        
        
        %Make each DFF per stim graph, one piece at a time
        %Note*: need to take off alpha value prior to exporting for AI
        %processing!
        %DFF
        figure(4*iStim-3);hold on;
        h = shadedErrorBar(time_vect,DFF_avg_split/max(DFF_avg_split),DFF_err_split/max(DFF_avg_split),...
            {'-','Color',colConst(iConst,:)});
        hold off;
        legHandDFF_ts = [legHandDFF_ts,h.mainLine];
        %SNR
        figure(4*iStim-2);hold on;
        h = shadedErrorBar(time_vect,SNR_avg_split/max(SNR_avg_split),SNR_err_split/max(SNR_avg_split),...
            {'-','Color',colConst(iConst,:)});
        hold off;
        legHandSNR_ts = [legHandSNR_ts,h.mainLine];
        
        %Get a vector of the maximum DFF for the histogram
        DFF_maxHold_split = max(DFF_mat_split(goodTimesTF,:),[],1);
        DFF_maxHist_split = nancat(1,DFF_maxHist_split,DFF_maxHold_split);
        %SNR
        SNR_maxHold_split = max(SNR_mat_split(goodTimesTF,:),[],1);
        SNR_maxHist_split = nancat(1,SNR_maxHist_split,SNR_maxHold_split);
        
        %Build a 3D vector containg the tHalfDecay Values
%         thisThalf = tHalfDecayCalc(DFF_mat_split,time_vect);
%         thisThalf = reshape(thisThalf,[1,1,length(thisThalf)]);
%         thisThalf = tHalf(
        
    end %iConst
    
    %add legend and labels to DFF time course broken out by construct
    figure(4*iStim - 3);
    legend(legHandDFF_ts,cellstr(constLev));
    xlabel('time (sec)');
    ylabel('\DeltaF/F');
    titleStr = sprintf(' %sFP',char(stimLev_split(iStim)));
    title(['\DeltaF/F values for',titleStr]);
    
    %add legend and labels to SNR time course broken out by construct
    figure(4*iStim - 2);
    legend(legHandSNR_ts,cellstr(constLev));
    xlabel('time (sec)');
    ylabel('SNR');
    titleStr = sprintf(' %sFP',char(stimLev_split(iStim)));
    title(['SNR values for',titleStr]);
    
    %Remove outliers for histogram data
    %DFF
    avgDFF = nanmean(DFF_maxHist_split,2);
    stdDFF = nanstd(DFF_maxHist_split,[],2);
    thold = repmat(avgDFF + 8*stdDFF,1,size(DFF_maxHist_split,2));
    DFF_maxHist_split(DFF_maxHist_split>thold) = NaN;
    %SNR
    avgSNR = nanmean(SNR_maxHist_split,2);
    stdSNR = nanstd(SNR_maxHist_split,[],2);
    thold = repmat(avgSNR + 8*stdSNR,1,size(SNR_maxHist_split,2));
    SNR_maxHist_split(SNR_maxHist_split>thold) = NaN;
    
 %%   %Building the scatterhist plots
    %Remove values of F0 higher than 9 from dataset as these always throw
    %off the fits
    badIdx = F0_mat_split>9;
    badIdx = logical(sum(badIdx,1));
    F0_mat_reduced = F0_mat_split;
    F0_mat_reduced(:,badIdx) = [];
    DFF_mat_reduced = DFF_maxHist_split;
    DFF_mat_reduced(:,badIdx) = [];
    SNR_mat_reduced = SNR_maxHist_split;
    SNR_mat_reduced(:,badIdx) = [];
    
    %Sparsify the number of dudes in the scatterhist plots. Fits and such
    %will be to these values as well
    skip = 30;
    %Sparsen the data points
    group = repmat(cellstr(constLev),size(F0_mat_reduced,2),1)';
    group_sub = group(:,1:skip:end);
    F0_mat_sub = F0_mat_reduced(:,1:skip:end);
    DFF_maxHist_sub = DFF_mat_reduced(:,1:skip:end);
    SNR_maxHist_sub = SNR_mat_reduced(:,1:skip:end);
    
    %Calculate R^2 stats for each dataset
    x_G6s = [F0_mat_sub(1,:)',ones(size(F0_mat_sub,2),1)];
    y_G6s_DFF = DFF_maxHist_sub(1,:)';
    y_G6s_SNR = SNR_maxHist_sub(1,:)';
    x_GAP = [F0_mat_sub(2,:)',ones(size(F0_mat_sub,2),1)];
    y_GAP_DFF = DFF_maxHist_sub(2,:)';
    y_GAP_SNR = SNR_maxHist_sub(2,:)';
    %Calculating stats
    [~,~,~,~,stat_G6s_DFF] = regress(y_G6s_DFF,x_G6s);
    [~,~,~,~,stat_G6s_SNR] = regress(y_G6s_SNR,x_G6s);
    [~,~,~,~,stat_GAP_DFF] = regress(y_GAP_DFF,x_GAP);
    [~,~,~,~,stat_GAP_SNR] = regress(y_GAP_SNR,x_GAP);
    
    %Calculating bin size to enforce bin size equality
    %F0
    F0_range = range(F0_mat_sub');
    xBinScale = 20;
    Nbin_F0 = [xBinScale,round(F0_range(2)/F0_range(1)*xBinScale)];
    %DFF
    yBinScale = 20;
    DFF_range = range(DFF_maxHist_sub');
    Nbin_DFF = [yBinScale,round(DFF_range(2)/DFF_range(1)*yBinScale)];
    %SNR
    SNR_range = range(SNR_maxHist_sub');
    Nbin_SNR = [yBinScale,round(SNR_range(2)/SNR_range(1)*yBinScale)];
    
    %Reshape variables for the scatterhist plots
    sizer = size(F0_mat_sub);
    F0_sc = reshape(F0_mat_sub,sizer(1)*sizer(2),1);
    DFF_sc = reshape(DFF_maxHist_sub,sizer(1)*sizer(2),1);
    SNR_sc = reshape(SNR_maxHist_sub,sizer(1)*sizer(2),1);
    group_sc = reshape(group_sub,sizer(1)*sizer(2),1);
    
    %Making the plots, start with DFF
    figure(4*iStim - 1);
    scatterhist(F0_sc,DFF_sc,...
        'Group',group_sc,...
        'Nbins',{Nbin_F0,Nbin_DFF},...
        'Style','Stairs');
    %legend(legHandDFF_ts,cellstr(constLev));
    xlabel('F0')
    ylabel('\DeltaF/F')
    titleStr = sprintf(' %sFP',char(stimLev_split(iStim)));
    title(['F0 v \DeltaF/F scatter for',titleStr]);
    %Add linear fits
    lsline(gca)
    %Add R^2 statistic
    text(6,1,['R^2 = ',num2str(stat_G6s_DFF(1))],'Color','b');
    text(6,1.5,['R^2 = ',num2str(stat_GAP_DFF(1))],'Color','r');
    
    %add legend and labels to F0 v SNR scattterplot
    figure(4*iStim);
    scatterhist(F0_sc,SNR_sc,...
        'Group',group_sc,...
        'Nbins',{Nbin_F0,Nbin_SNR},...
        'Style','Stairs');
    %legend(legHandSNR_ts,cellstr(constLev));
    xlabel('F0')
    ylabel('SNR')
    titleStr = sprintf(' %sFP',char(stimLev_split(iStim)));
    title(['F0 v SNR scatter for',titleStr]);
    %Add linear fits
    lsline(gca)
    %Add R^2 statistic
    text(6,10,['R^2 = ',num2str(stat_G6s_SNR(1))],'Color','b');
    text(6,20,['R^2 = ',num2str(stat_GAP_SNR(1))],'Color','r');
    
%%    %make figure for the comparison of histograms
    %make figure DFF
    figure(2*iStim + 199);hold on;
    for iHist = 1:numel(constLev)
        histogram(DFF_maxHist_split(iHist,:),100,...
            'Normalization','pdf',...
            'DisplayStyle','Stairs');
        xlabel('DFF')
        ylabel('Binned ROI counts')
        legend(cellstr(constLev));
        title(['\DeltaF/F histogram for',titleStr]);
        colormap(hsv(numel(constLev)))
    end
    %add median lines DFF
    plotYlim = get(gca,'Ylim');
    plot([nanmean(DFF_maxHist_split(1,:)),nanmean(DFF_maxHist_split(1,:))],...
        [plotYlim(1),plotYlim(2)]);
    dffStr = sprintf('G6s = %0.2f',nanmean(DFF_maxHist_split(1,:)));
    text(nanmean(DFF_maxHist_split(1,:)),4*(plotYlim(1)+plotYlim(2))/5,dffStr);
    plot([nanmean(DFF_maxHist_split(2,:)),nanmean(DFF_maxHist_split(2,:))],...
        [plotYlim(1),plotYlim(2)]);
    dffStr = sprintf('GAP43-G6s = %0.2f',nanmean(DFF_maxHist_split(2,:)));
    text(nanmean(DFF_maxHist_split(2,:)),2*(plotYlim(1)+plotYlim(2))/3,dffStr);
    %make figure SNR
    figure(2*iStim + 200);hold on;
    for iHist = 1:numel(constLev)
        histogram(SNR_maxHist_split(iHist,:),100,...
            'Normalization','pdf',...
            'DisplayStyle','Stairs');
        xlabel('SNR')
        ylabel('Binned ROI counts')
        legend(cellstr(constLev));
        title(['SNR histogram for',titleStr]);
        colormap(hsv(numel(constLev)))
    end
    %add median lines SNR
    plotYlim = get(gca,'Ylim');
    plot([nanmean(SNR_maxHist_split(1,:)),nanmean(SNR_maxHist_split(1,:))],...
        [plotYlim(1),plotYlim(2)]);
    SNRStr = sprintf('G6s = %0.2f',nanmean(SNR_maxHist_split(1,:)));
    text(nanmean(SNR_maxHist_split(1,:)),4*(plotYlim(1)+plotYlim(2))/5,SNRStr);
    plot([nanmean(SNR_maxHist_split(2,:)),nanmean(SNR_maxHist_split(2,:))],...
        [plotYlim(1),plotYlim(2)]);
    SNRStr = sprintf('GAP43-G6s = %0.2f',nanmean(SNR_maxHist_split(2,:)));
    text(nanmean(SNR_maxHist_split(2,:)),2*(plotYlim(1)+plotYlim(2))/3,SNRStr);
    
end %iStim
%% make full summary figures

%have to convert to string mat, then cell array from nominal variables for
%legends
%make laser plot
figure;
cmap = lines(numel(constLev));
for i = 1:numel(constLev)
    thisRedF_split = RedF_split(:,:,i);
    hold on;
    h = histogram(RedF_split(:,:,i),50,...
    'Normalization','pdf',...k
    'DisplayStyle','Stairs');
    h.EdgeColor = cmap(i,:);
    %Get the median value for each distriibtution, add to graph
    thisMed = nanmedian(thisRedF_split(:));
    ylims = ylim;
    plot([thisMed,thisMed],ylims,'Color',cmap(i,:));
    thisMedText = num2str(thisMed);
    text(thisMed,4/(i*5)*(ylims(2) - ylims(1)),thisMedText,'Color',cmap(i,:));
    %Labels
    xlabel('Red fluorescence intensity (A.U.)')
    ylabel('PDF')
    title('Distribution of Red Fluorescence per ROI for each construct')
    legend(cellstr(constLev));
    
end

%make ROI sizes plot
figure;
cmap = lines(numel(constLev));
for i = 1:numel(constLev)
    thisROIsizes_split = ROIsizes_split(:,:,i);
    hold on;
    h = histogram(ROIsizes_split(:,:,i),50,...
    'Normalization','pdf',...
    'DisplayStyle','Stairs');
    h.EdgeColor = cmap(i,:);
    %Get the median value for each distriibtution, add to graph
    thisMed = nanmedian(thisROIsizes_split(:));
    ylims = ylim;
    plot([thisMed,thisMed],ylims,'Color',cmap(i,:));
    thisMedText = num2str(thisMed);
    text(thisMed,4/(i*5)*(ylims(2) - ylims(1)),thisMedText,'Color',cmap(i,:));
    %Labels
    xlabel('Size of ROIs (\mum^2)')
    ylabel('Binned ROI counts')
    title('Distribution of ROI sizes for each construct')
    legend(cellstr(constLev));
    
end

%make DFF summary plot
set(0,'DefaultAxesColorOrder',lines(numel(constLev)))
figure;
stimChar = char(stimLev);
str = [];for i = 1:size(stimLev,2);str = [str,' ',stimChar(i,:),'FP'];end
xplot = repmat(sscanf(str,['%d%*s']),[1,size(DFF_max_split,2)]);
errorbar(xplot,DFF_max_all,DFF_errMat_all);
legend(cellstr(constLev),'Location','NorthWest');
%text(xplot(3:end,numel(constLev)),DFF_max(3:end,numel(constLev))+1,'*')
xlabel('Field Potentials')
ylabel('\DeltaF/F')
colormap(lines(numel(constLev)))
title('\DeltaF/F population summary (total mask, error bar sem)')

%make SNR summary plot
set(0,'DefaultAxesColorOrder',lines(numel(constLev)))
figure;
stimChar = char(stimLev);
str = [];for i = 1:size(stimLev,2);str = [str,' ',stimChar(i,:),'FP'];end
xplot = repmat(sscanf(str,['%d%*s']),[1,size(DFF_max_split,2)]);
errorbar(xplot,SNR_max_split,DFF_errMat_split);
legend(cellstr(constLev),'Location','NorthWest');
%text(xplot(3:end,numel(constLev)),DFF_max(3:end,numel(constLev))+1,'*')
xlabel('Field Potentials')
ylabel('\DeltaF/F')
colormap(lines(numel(constLev)))
title('SNR population summary (split masks, error bars std)')